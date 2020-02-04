"""
Utility functions for dealing with input args to the sourmash command line.
"""
import sys
import os
import argparse
from . import signature
from .logging import notify, error

from .index import LinearIndex
from . import signature as sig
from .sbt import SBT
from .sbtmh import SigLeaf
from .lca import lca_utils
import sourmash

DEFAULT_LOAD_K = 31
DEFAULT_N = 500


def get_moltype(sig, require=False):
    if sig.minhash.is_molecule_type('DNA'):
        moltype = 'DNA'
    elif sig.minhash.is_molecule_type('dayhoff'):
        moltype = 'dayhoff'
    elif sig.minhash.is_molecule_type('hp'):
        moltype = 'hp'
    elif sig.minhash.is_molecule_type('protein'):
        moltype = 'protein'
    else:
        raise ValueError('unknown molecule type for sig {}'.format(sig.name()))

    return moltype


def get_sigtype(sig):
    if sig.minhash.scaled:
        return 'scaled'
    elif sig.minhash.num:
        return 'bounded'
    else:
        raise Exception("unknown signature type")


class CheckCompatibleSignatures(object):
    "Track compatibility of a bunch of signatures."
    def __init__(self):
        self.ksize = None
        self.moltype = None
        self.sigtype = None
        self.problem = None

    def check(self, new_sig):
        """Check new signature against previous ksize/moltype/sigtype.

        If this is the first signature, record it for future checks."""

        ksize = new_sig.minhash.ksize
        moltype = get_moltype(new_sig)
        sigtype = get_sigtype(new_sig)

        if self.ksize is None:
            self.ksize = ksize
        elif ksize != self.ksize:
            self.problem = 'ksize'
            raise ValueError("new ksize {}, vs previous {} ".format(ksize,
                                                                self.ksize))
        else:
            assert ksize == self.ksize

        if self.moltype is None:
            self.moltype = moltype
            self.problem = 'moltype'
        elif moltype != self.moltype:
            raise ValueError("new moltype {}, vs previous {}".format(moltype,
                                                                self.moltype))
        else:
            assert moltype == self.moltype

        if self.sigtype is None:
            self.sigtype = sigtype
        elif sigtype != self.sigtype:
            self.problem = 'sigtype'
            raise ValueError("new sigtype {}, vs previous {}".format(sigtype,
                                                                self.sigtype))
        else:
            assert sigtype == self.sigtype

    def check_list(self, siglist):
        "Check multiple signatures."
        for s in siglist:
            self.check(s)


def calculate_moltype(args, default=None):
    if args.protein:
        if args.dna is True:
            error('cannot specify both --dna/--rna and --protein!')
            sys.exit(-1)
        args.dna = False

    moltype = default
    if args.dna:
        moltype = 'DNA'
    elif args.dayhoff:
        moltype = 'dayhoff'
    elif args.hp:
        moltype = 'hp'
    elif args.protein:
        moltype = 'protein'

    return moltype


def load_query_signature(filename, ksize, select_moltype):
    try:
        sl = signature.load_signatures(filename,
                                       ksize=ksize,
                                       select_moltype=select_moltype,
                                       do_raise=True)
        sl = list(sl)
    except (IOError, ValueError):
        error("Cannot open file '{}'", filename)
        sys.exit(-1)

    if len(sl) and ksize is None:
        ksizes = set([ ss.minhash.ksize for ss in sl ])
        if len(ksizes) == 1:
            ksize = ksizes.pop()
            sl = [ ss for ss in sl if ss.minhash.ksize == ksize ]
            notify('select query k={} automatically.', ksize)
        elif DEFAULT_LOAD_K in ksizes:
            sl = [ ss for ss in sl if ss.minhash.ksize == DEFAULT_LOAD_K ]
            notify('selecting default query k={}.', DEFAULT_LOAD_K)
    elif ksize:
        notify('selecting specified query k={}', ksize)

    if len(sl) != 1:
        error('When loading query from "{}"', filename)
        error('{} signatures matching ksize and molecule type;', len(sl))
        error('need exactly one. Specify --ksize or --dna, --rna, or --protein.')
        sys.exit(-1)

    return sl[0]


class LoadSingleSignatures(object):
    def __init__(self, filelist,  ksize=None, select_moltype=None,
                 ignore_files=set()):
        self.filelist = filelist
        self.ksize = ksize
        self.select_moltype = select_moltype
        self.ignore_files = ignore_files

        self.skipped_ignore = 0
        self.skipped_nosig = 0
        self.ksizes = set()
        self.moltypes = set()

    def __iter__(self):
        for filename in self.filelist:
            if filename in self.ignore_files:
                self.skipped_ignore += 1
                continue

            sl = signature.load_signatures(filename,
                                           ksize=self.ksize,
                                           select_moltype=self.select_moltype)
            sl = list(sl)
            if len(sl) == 0:
                self.skipped_nosig += 1
                continue

            for query in sl:
                query_moltype = get_moltype(query)
                query_ksize = query.minhash.ksize

                self.ksizes.add(query_ksize)
                self.moltypes.add(query_moltype)

            if len(self.ksizes) > 1 or len(self.moltypes) > 1:
                raise ValueError('multiple k-mer sizes/molecule types present')

            for query in sl:
                yield filename, query, query_moltype, query_ksize


def traverse_find_sigs(dirnames, yield_all_files=False):
    for dirname in dirnames:
        if (dirname.endswith('.sig') or yield_all_files) and os.path.isfile(dirname):
            yield dirname
            continue

        for root, dirs, files in os.walk(dirname):
            for name in files:
                if name.endswith('.sig') or yield_all_files:
                    fullname = os.path.join(root, name)
                    yield fullname


def filter_compatible_signatures(query, siglist, force=False):
    for ss in siglist:
        if check_signatures_are_compatible(query, ss):
            yield ss
        else:
            if not force:
                raise ValueError("incompatible signature")


def check_signatures_are_compatible(query, subject):
    # is one scaled, and the other not? cannot do search
    if query.minhash.scaled and not subject.minhash.scaled or \
       not query.minhash.scaled and subject.minhash.scaled:
       error("signature {} and {} are incompatible - cannot compare.",
             query.name(), subject.name())
       if query.minhash.scaled:
           error("{} was calculated with --scaled, {} was not.",
                 query.name(), subject.name())
       if subject.minhash.scaled:
           error("{} was calculated with --scaled, {} was not.",
                 subject.name(), query.name())
       return 0

    return 1


def check_tree_is_compatible(treename, tree, query, is_similarity_query):
    # get a minhash from the tree
    leaf = next(iter(tree.leaves()))
    tree_mh = leaf.data.minhash

    query_mh = query.minhash

    if tree_mh.ksize != query_mh.ksize:
        error("ksize on tree '{}' is {};", treename, tree_mh.ksize)
        error('this is different from query ksize of {}.', query_mh.ksize)
        return 0

    # is one scaled, and the other not? cannot do search.
    if (tree_mh.scaled and not query_mh.scaled) or \
       (query_mh.scaled and not tree_mh.scaled):
        error("for tree '{}', tree and query are incompatible for search.",
              treename)
        if tree_mh.scaled:
            error("tree was calculated with scaled, query was not.")
        else:
            error("query was calculated with scaled, tree was not.")
        return 0

    # are the scaled values incompatible? cannot downsample tree for similarity
    if tree_mh.scaled and tree_mh.scaled < query_mh.scaled and \
      is_similarity_query:
        error("for tree '{}', scaled value is smaller than query.", treename)
        error("tree scaled: {}; query scaled: {}. Cannot do similarity search.",
              tree_mh.scaled, query_mh.scaled)
        return 0

    return 1


def load_dbs_and_sigs(filenames, query, is_similarity_query, traverse=False):
    """
    Load one or more SBTs, LCAs, and/or signatures.

    Check for compatibility with query.
    """
    query_ksize = query.minhash.ksize
    query_moltype = get_moltype(query)

    n_signatures = 0
    n_databases = 0
    databases = []
    for sbt_or_sigfile in filenames:
        notify('loading from {}...', sbt_or_sigfile, end='\r')
        # are we collecting signatures from a directory/path?
        if traverse and os.path.isdir(sbt_or_sigfile):
            for sigfile in traverse_find_sigs([sbt_or_sigfile]):
                try:
                    siglist = sig.load_signatures(sigfile,
                                                  ksize=query_ksize,
                                                  select_moltype=query_moltype)
                    siglist = filter_compatible_signatures(query, siglist, 1)
                    linear = LinearIndex(siglist, filename=sigfile)
                    databases.append((linear, sbt_or_sigfile, False))
                    notify('loaded {} signatures from {}', len(linear),
                           sigfile, end='\r')
                    n_signatures += len(linear)
                except Exception:                       # ignore errors with traverse
                    pass

            # done! jump to beginning of main 'for' loop
            continue

        # no traverse? try loading as an SBT.
        try:
            tree = SBT.load(sbt_or_sigfile, leaf_loader=SigLeaf.load)

            if not check_tree_is_compatible(sbt_or_sigfile, tree, query,
                                            is_similarity_query):
                sys.exit(-1)

            databases.append((tree, sbt_or_sigfile, 'SBT'))
            notify('loaded SBT {}', sbt_or_sigfile, end='\r')
            n_databases += 1

            # done! jump to beginning of main 'for' loop
            continue
        except (ValueError, EnvironmentError):
            # not an SBT - try as an LCA
            pass

        # ok. try loading as an LCA.
        try:
            lca_db = lca_utils.LCA_Database()
            lca_db.load(sbt_or_sigfile)

            assert query_ksize == lca_db.ksize
            query_scaled = query.minhash.scaled

            notify('loaded LCA {}', sbt_or_sigfile, end='\r')
            n_databases += 1

            databases.append((lca_db, sbt_or_sigfile, 'LCA'))

            continue
        except (ValueError, TypeError, EnvironmentError):
            # not an LCA database - try as a .sig
            pass

        # not a tree? try loading as a signature.
        try:
            siglist = sig.load_signatures(sbt_or_sigfile,
                                          ksize=query_ksize,
                                          select_moltype=query_moltype)
            siglist = list(siglist)
            if len(siglist) == 0:         # file not found, or parse error?
                raise ValueError

            siglist = filter_compatible_signatures(query, siglist, False)
            linear = LinearIndex(siglist, filename=sbt_or_sigfile)
            databases.append((linear, sbt_or_sigfile, 'signature'))

            notify('loaded {} signatures from {}', len(linear),
                   sbt_or_sigfile, end='\r')
            n_signatures += len(linear)
        except (EnvironmentError, ValueError):
            error("\nCannot open file '{}'", sbt_or_sigfile)
            sys.exit(-1)

    notify(' '*79, end='\r')
    if n_signatures and n_databases:
        notify('loaded {} signatures and {} databases total.', n_signatures, 
                                                               n_databases)
    elif n_signatures:
        notify('loaded {} signatures.', n_signatures)
    elif n_databases:
        notify('loaded {} databases.', n_databases)
    else:
        sys.exit(-1)

    if databases:
        print('')

    return databases


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
    def __init__(self, filename, mode='wt'):
        self.filename = filename
        self.mode = mode
        self.fp = None

    def open(self):
        if self.filename == '-' or self.filename is None:
            return sys.stdout
        self.fp = open(self.filename, self.mode)
        return self.fp

    def __enter__(self):
        return self.open()

    def __exit__(self, type, value, traceback):
        # do we need to handle exceptions here?
        if self.fp:
            self.fp.close()

        return False
