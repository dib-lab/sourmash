"""
Utility functions for sourmash CLI commands.
"""
import sys
import os
import argparse
import itertools
from enum import Enum
import traceback
import gzip
import zipfile

import screed

from sourmash.sbtmh import load_sbt_index
from sourmash.lca.lca_db import load_single_database
import sourmash.exceptions

from . import signature
from .logging import notify, error, debug_literal

from .index import (LinearIndex, ZipFileLinearIndex, MultiIndex)
from . import signature as sig
from .sbt import SBT
from .sbtmh import SigLeaf
from .lca import LCA_Database
import sourmash

DEFAULT_LOAD_K = 31


def get_moltype(sig, require=False):
    mh = sig.minhash
    if mh.moltype in ('DNA', 'dayhoff', 'hp', 'protein'):
        moltype = mh.moltype
    else:
        raise ValueError('unknown molecule type for sig {}'.format(sig))

    return moltype


def calculate_moltype(args, default=None):
    moltype = default

    n = 0
    if args.dna:
        moltype = 'DNA'
        n += 1
    if args.dayhoff:
        moltype = 'dayhoff'
        n += 1
    if args.hp:
        moltype = 'hp'
        n += 1
    if args.protein:
        moltype = 'protein'
        n += 1

    if n > 1:
        error("cannot specify more than one of --dna/--rna/--protein/--hp/--dayhoff")
        sys.exit(-1)

    return moltype


def load_query_signature(filename, ksize, select_moltype, select_md5=None):
    """Load a single signature to use as a query.

    Uses load_file_as_signatures underneath, so can load from collections
    and indexed databases.
    """
    try:
        sl = load_file_as_signatures(filename, ksize=ksize,
                                     select_moltype=select_moltype)
        sl = list(sl)
    except (OSError, ValueError):
        error(f"Cannot open query file '{filename}'")
        sys.exit(-1)

    if len(sl) and select_md5:
        found_sig = None
        for sig in sl:
            sig_md5 = sig.md5sum()
            if sig_md5.startswith(select_md5.lower()):
                # make sure we pick only one --
                if found_sig is not None:
                    error(f"Error! Multiple signatures start with md5 '{select_md5}'")
                    error("Please use a longer --md5 selector.")
                    sys.exit(-1)
                else:
                    found_sig = sig

            sl = [found_sig]

    if len(sl) and ksize is None:
        ksizes = set([ ss.minhash.ksize for ss in sl ])
        if len(ksizes) == 1:
            ksize = ksizes.pop()
            sl = [ ss for ss in sl if ss.minhash.ksize == ksize ]
            notify(f'select query k={ksize} automatically.')
        elif DEFAULT_LOAD_K in ksizes:
            sl = [ ss for ss in sl if ss.minhash.ksize == DEFAULT_LOAD_K ]
            notify(f'selecting default query k={DEFAULT_LOAD_K}.')
    elif ksize:
        notify(f'selecting specified query k={ksize}')

    if len(sl) != 1:
        error(f"When loading query from '{filename}'", filename)
        error(f'{len(sl)} signatures matching ksize and molecule type;')
        error('need exactly one. Specify --ksize or --dna, --rna, or --protein.')
        sys.exit(-1)

    return sl[0]


def _check_suffix(filename, endings):
    for ending in endings:
        if filename.endswith(ending):
            return True
    return False


def traverse_find_sigs(filenames, yield_all_files=False):
    """Find all .sig and .sig.gz files in & beneath 'filenames'.

    By default, this function returns files with .sig and .sig.gz extensions.
    If 'yield_all_files' is True, this will return _all_ files
    (but not directories).
    """
    endings = ('.sig', '.sig.gz')
    for filename in filenames:
        # check for files in filenames:
        if os.path.isfile(filename):
            if yield_all_files or _check_suffix(filename, endings):
                yield filename

        # filename is a directory -- traverse beneath!
        elif os.path.isdir(filename):
            for root, dirs, files in os.walk(filename):
                for name in files:
                    fullname = os.path.join(root, name)
                    if yield_all_files or _check_suffix(fullname, endings):
                        yield fullname


def load_dbs_and_sigs(filenames, query, is_similarity_query, *, cache_size=None):
    """
    Load one or more SBTs, LCAs, and/or collections of signatures.

    Check for compatibility with query.

    This is basically a user-focused wrapping of _load_databases.
    """
    query_mh = query.minhash

    containment = True
    if is_similarity_query:
        containment = False

    databases = []
    for filename in filenames:
        notify(f'loading from {filename}...', end='\r')

        try:
            db = _load_database(filename, False, cache_size=cache_size)
        except ValueError as e:
            # cannot load database!
            notify(str(e))
            sys.exit(-1)

        try:
            db = db.select(moltype=query_mh.moltype,
                           ksize=query_mh.ksize,
                           num=query_mh.num,
                           scaled=query_mh.scaled,
                           containment=containment)
        except ValueError as exc:
            # incompatible collection specified!
            notify(f"ERROR: cannot use '{filename}' for this query.")
            notify(str(exc))
            sys.exit(-1)

        # 'select' returns nothing => all signatures filtered out. fail!
        if not db:
            notify(f"no compatible signatures found in '{filename}'")
            sys.exit(-1)

        databases.append(db)

    # calc num loaded info.
    n_signatures = 0
    n_databases = 0
    for db in databases:
        if db.is_database:
            n_databases += 1
        else:
            n_signatures += len(db)

    notify(' '*79, end='\r')
    if n_signatures and n_databases:
        notify(f'loaded {n_signatures} signatures and {n_databases} databases total.')
    elif n_signatures and not n_databases:
        notify(f'loaded {n_signatures} signatures.')
    elif n_databases and not n_signatures:
        notify(f'loaded {n_databases} databases.')

    if databases:
        print('')
    else:
        notify('** ERROR: no signatures or databases loaded?')
        sys.exit(-1)

    return databases


def _load_stdin(filename, **kwargs):
    "Load collection from .sig file streamed in via stdin"
    db = None
    if filename == '-':
        db = LinearIndex.load(sys.stdin)

    return db


def _multiindex_load_from_pathlist(filename, **kwargs):
    "Load collection from a list of signature/database files"
    db = MultiIndex.load_from_pathlist(filename)

    return db


def _multiindex_load_from_path(filename, **kwargs):
    "Load collection from a directory."
    traverse_yield_all = kwargs['traverse_yield_all']
    db = MultiIndex.load_from_path(filename, traverse_yield_all)

    return db


def _load_sigfile(filename, **kwargs):
    "Load collection from a signature JSON file"
    try:
        db = LinearIndex.load(filename)
    except sourmash.exceptions.SourmashError as exc:
        raise ValueError(exc)

    return db


def _load_sbt(filename, **kwargs):
    "Load collection from an SBT."
    cache_size = kwargs.get('cache_size')

    try:
        db = load_sbt_index(filename, cache_size=cache_size)
    except (FileNotFoundError, TypeError) as exc:
        raise ValueError(exc)

    return db


def _load_revindex(filename, **kwargs):
    "Load collection from an LCA database/reverse index."
    db, _, _ = load_single_database(filename)
    return db


def _load_zipfile(filename, **kwargs):
    "Load collection from a .zip file."
    db = None
    if filename.endswith('.zip'):
        traverse_yield_all = kwargs['traverse_yield_all']
        db = ZipFileLinearIndex.load(filename,
                                     traverse_yield_all=traverse_yield_all)
    return db


# all loader functions, in order.
_loader_functions = [
    ("load from stdin", _load_stdin),
    ("load from directory", _multiindex_load_from_path),
    ("load from sig file", _load_sigfile),
    ("load from file list", _multiindex_load_from_pathlist),
    ("load SBT", _load_sbt),
    ("load revindex", _load_revindex),
    ("load collection from zipfile", _load_zipfile),
    ]


def _load_database(filename, traverse_yield_all, *, cache_size=None):
    """Load file as a database - list of signatures, LCA, SBT, etc.

    Return Index object.

    This is an internal function used by other functions in sourmash_args.
    """
    loaded = False

    # iterate through loader functions, trying them all. Catch ValueError
    # but nothing else.
    for (desc, load_fn) in _loader_functions:
        try:
            debug_literal(f"_load_databases: trying loader fn {desc}")
            db = load_fn(filename,
                         traverse_yield_all=traverse_yield_all,
                         cache_size=cache_size)
        except ValueError as exc:
            debug_literal(f"_load_databases: FAIL on fn {desc}.")
            debug_literal(traceback.format_exc())

        if db:
            loaded = True
            break

    # check to see if it's a FASTA/FASTQ record (i.e. screed loadable)
    # so we can provide a better error message to users.
    if not loaded:
        successful_screed_load = False
        it = None
        try:
            # CTB: could be kind of time consuming for a big record, but at the
            # moment screed doesn't expose format detection cleanly.
            with screed.open(filename) as it:
                record = next(iter(it))
            successful_screed_load = True
        except:
            pass

        if successful_screed_load:
            raise ValueError(f"Error while reading signatures from '{filename}' - got sequences instead! Is this a FASTA/FASTQ file?")

    if not loaded:
        raise ValueError(f"Error while reading signatures from '{filename}'.")

    if loaded:                  # this is a bit redundant but safe > sorry
        assert db

    return db


def load_file_as_index(filename, yield_all_files=False):
    """Load 'filename' as a database; generic database loader.

    If 'filename' contains an SBT or LCA indexed database, or a regular
    Zip file, will return the appropriate objects. If a Zip file and
    yield_all_files=True, will try to load all files within zip, not just
    .sig files.

    If 'filename' is a JSON file containing one or more signatures, will
    return an Index object containing those signatures.

    If 'filename' is a directory, will load *.sig underneath
    this directory into an Index object. If yield_all_files=True, will
    attempt to load all files.
    """
    return _load_database(filename, yield_all_files)


def load_file_as_signatures(filename, select_moltype=None, ksize=None,
                            yield_all_files=False,
                            progress=None):
    """Load 'filename' as a collection of signatures. Return an iterable.

    If 'filename' contains an SBT or LCA indexed database, or a regular
    Zip file, will return a signatures() generator. If a Zip file and
    yield_all_files=True, will try to load all files within zip, not just
    .sig files.

    If 'filename' is a JSON file containing one or more signatures, will
    return a list of those signatures.

    If 'filename' is a directory, will load *.sig
    underneath this directory into a list of signatures. If
    yield_all_files=True, will attempt to load all files.

    Applies selector function if select_moltype and/or ksize are given.
    """
    if progress:
        progress.notify(filename)

    db = _load_database(filename, yield_all_files)
    db = db.select(moltype=select_moltype, ksize=ksize)
    loader = db.signatures()

    if progress:
        return progress.start_file(filename, loader)
    else:
        return loader


def load_pathlist_from_file(filename):
    "Load a list-of-files text file."
    try:
        with open(filename, 'rt') as fp:
            file_list = [ x.rstrip('\r\n') for x in fp ]
        file_list = set(file_list)
        if not file_list:
            raise ValueError("pathlist is empty")
        for checkfile in file_list:
            if not os.path.exists(checkfile):
                raise ValueError(f"file '{checkfile}' listed inside the pathlist contains a badly formatted file")
    except IOError:
        raise ValueError(f"pathlist file '{filename}' does not exist")    
    except OSError:
        raise ValueError(f"cannot open file '{filename}'")
    except UnicodeDecodeError:
        raise ValueError(f"cannot parse file '{filename}' as list of filenames")
    return file_list


class FileOutput(object):
    """A context manager for file outputs that handles sys.stdout gracefully.

    Usage:

       with FileOutput(filename, mode) as fp:
          ...

    does what you'd expect, but it handles the situation where 'filename'
    is '-' or None. This makes it nicely compatible with argparse usage,
    e.g.

    p = argparse.ArgumentParser()
    p.add_argument('--output')
    args = p.parse_args()
    ...
    with FileOutput(args.output, 'wt') as fp:
       ...

    will properly handle no argument or '-' as sys.stdout.
    """
    def __init__(self, filename, mode='wt', *, newline=None, encoding='utf-8'):
        self.filename = filename
        self.mode = mode
        self.fp = None
        self.newline = newline
        self.encoding = encoding

    def open(self):
        if self.filename == '-' or self.filename is None:
            return sys.stdout
        self.fp = open(self.filename, self.mode, newline=self.newline,
                       encoding=self.encoding)
        return self.fp

    def __enter__(self):
        return self.open()

    def __exit__(self, type, value, traceback):
        # do we need to handle exceptions here?
        if self.fp:
            self.fp.close()

        return False

class FileOutputCSV(FileOutput):
    """A context manager for CSV file outputs.

    Usage:

       with FileOutputCSV(filename) as fp:
          ...

    does what you'd expect, but it handles the situation where 'filename'
    is '-' or None. This makes it nicely compatible with argparse usage,
    e.g.

    p = argparse.ArgumentParser()
    p.add_argument('--output')
    args = p.parse_args()
    ...
    with FileOutputCSV(args.output) as w:
       ...

    will properly handle no argument or '-' as sys.stdout.
    """
    def __init__(self, filename):
        self.filename = filename
        self.fp = None

    def open(self):
        if self.filename == '-' or self.filename is None:
            return sys.stdout
        self.fp = open(self.filename, 'w', newline='')
        return self.fp


class SignatureLoadingProgress(object):
    """A wrapper for signature loading progress reporting.

    Instantiate this class once, and then pass it to load_file_as_signatures
    with progress=<obj>.

    Alternatively, call obj.start_file(filename, iter) each time you
    start loading signatures from a new file via iter.

    You can optionally notify of reading a file with `.notify(filename)`.
    """
    def __init__(self, reporting_interval=10):
        self.n_sig = 0
        self.interval = reporting_interval
        self.screen_width = 79

    def short_notify(self, msg_template, *args, **kwargs):
        """Shorten the notification message so that it fits on one line.

        Good for repeating notifications with end='\r' especially...
        """

        msg = msg_template.format(*args, **kwargs)
        end = kwargs.get('end', '\n')
        w = self.screen_width

        if len(msg) > w:
            truncate_len = len(msg) - w + 3
            msg = '<<<' + msg[truncate_len:]

        notify(msg, end=end)

    def notify(self, filename):
        self.short_notify("...reading from file '{}'",
                          filename, end='\r')

    def start_file(self, filename, loader):
        n_this = 0
        n_before = self.n_sig

        try:
            for result in loader:
                # track n from this file, as well as total n
                n_this += 1
                n_total = n_before + n_this
                if n_this and n_total % self.interval == 0:
                    self.short_notify("...loading from '{}' / {} sigs total",
                                      filename, n_total, end='\r')

                yield result
        except KeyboardInterrupt:
            # might as well nicely handle CTRL-C while we're at it!
            notify('\n(CTRL-C received! quitting.)')
            sys.exit(-1)
        finally:
            self.n_sig += n_this

        self.short_notify("loaded {} sigs from '{}'", n_this, filename)


#
# enum and classes for saving signatures progressively
#

class _BaseSaveSignaturesToLocation:
    "Base signature saving class. Track location (if any) and count."
    def __init__(self, location):
        self.location = location
        self.count = 0

    def __repr__(self):
        raise NotImplementedError

    def __len__(self):
        return self.count

    def __enter__(self):
        "provide context manager functionality"
        self.open()
        return self

    def __exit__(self, type, value, traceback):
        "provide context manager functionality"
        self.close()

    def add(self, ss):
        self.count += 1


class SaveSignatures_NoOutput(_BaseSaveSignaturesToLocation):
    "Do not save signatures."
    def __repr__(self):
        return 'SaveSignatures_NoOutput()'

    def open(self):
        pass

    def close(self):
        pass


class SaveSignatures_Directory(_BaseSaveSignaturesToLocation):
    "Save signatures within a directory, using md5sum names."
    def __init__(self, location):
        super().__init__(location)
        
    def __repr__(self):
        return f"SaveSignatures_Directory('{self.location}')"

    def close(self):
        pass

    def open(self):
        try:
            os.mkdir(self.location)
        except FileExistsError:
            pass
        except:
            notify("ERROR: cannot create signature output directory '{}'",
                   self.location)
            sys.exit(-1)

    def add(self, ss):
        super().add(ss)
        md5 = ss.md5sum()

        # don't overwrite even if duplicate md5sum
        outname = os.path.join(self.location, f"{md5}.sig.gz")
        if os.path.exists(outname):
            i = 0
            while 1:
                outname = os.path.join(self.location, f"{md5}_{i}.sig.gz")
                if not os.path.exists(outname):
                    break
                i += 1

        with gzip.open(outname, "wb") as fp:
            sig.save_signatures([ss], fp, compression=1)


class SaveSignatures_SigFile(_BaseSaveSignaturesToLocation):
    "Save signatures within a directory, using md5sum names."
    def __init__(self, location):
        super().__init__(location)
        self.keep = []
        self.compress = 0
        if self.location.endswith('.gz'):
            self.compress = 1

    def __repr__(self):
        return f"SaveSignatures_SigFile('{self.location}')"

    def open(self):
        pass

    def close(self):
        if self.location == '-':
            sourmash.save_signatures(self.keep, sys.stdout)
        else:
            # text mode? encode in utf-8
            mode = "w"
            encoding = 'utf-8'

            # compressed? bytes & binary.
            if self.compress:
                encoding = None
                mode = "wb"

            with open(self.location, mode, encoding=encoding) as fp:
                sourmash.save_signatures(self.keep, fp,
                                         compression=self.compress)

    def add(self, ss):
        super().add(ss)
        self.keep.append(ss)


class SaveSignatures_ZipFile(_BaseSaveSignaturesToLocation):
    "Save compressed signatures in an uncompressed Zip file."
    def __init__(self, location):
        super().__init__(location)
        self.zf = None
        
    def __repr__(self):
        return f"SaveSignatures_ZipFile('{self.location}')"

    def close(self):
        self.zf.close()

    def open(self):
        self.zf = zipfile.ZipFile(self.location, 'w', zipfile.ZIP_STORED)

    def _exists(self, name):
        try:
            self.zf.getinfo(name)
            return True
        except KeyError:
            return False

    def add(self, ss):
        assert self.zf
        super().add(ss)

        md5 = ss.md5sum()
        outname = f"signatures/{md5}.sig.gz"

        # don't overwrite even if duplicate md5sum.
        if self._exists(outname):
            i = 0
            while 1:
                outname = os.path.join(self.location, f"{md5}_{i}.sig.gz")
                if not self._exists(outname):
                    break
                i += 1

        json_str = sourmash.save_signatures([ss], compression=1)
        self.zf.writestr(outname, json_str)


class SigFileSaveType(Enum):
    SIGFILE = 1
    SIGFILE_GZ = 2
    DIRECTORY = 3
    ZIPFILE = 4
    NO_OUTPUT = 5

_save_classes = {
    SigFileSaveType.SIGFILE: SaveSignatures_SigFile,
    SigFileSaveType.SIGFILE_GZ: SaveSignatures_SigFile,
    SigFileSaveType.DIRECTORY: SaveSignatures_Directory,
    SigFileSaveType.ZIPFILE: SaveSignatures_ZipFile,
    SigFileSaveType.NO_OUTPUT: SaveSignatures_NoOutput
}


def SaveSignaturesToLocation(filename, *, force_type=None):
    """Create and return an appropriate object for progressive saving of
    signatures."""
    save_type = None
    if not force_type:
        if filename is None:
            save_type = SigFileSaveType.NO_OUTPUT
        elif filename.endswith('/'):
            save_type = SigFileSaveType.DIRECTORY
        elif filename.endswith('.gz'):
            save_type = SigFileSaveType.SIGFILE_GZ
        elif filename.endswith('.zip'):
            save_type = SigFileSaveType.ZIPFILE
        else:
            # default to SIGFILE intentionally!
            save_type = SigFileSaveType.SIGFILE
    else:
        save_type = force_type

    cls = _save_classes.get(save_type)
    if cls is None:
        raise Exception("invalid save type; this should never happen!?")

    return cls(filename)
