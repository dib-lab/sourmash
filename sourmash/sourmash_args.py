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


class SignatureParams(object):
    def __init__(self, ksizes, moltypes, num_vals, scaled_vals):
        self.ksizes = set(ksizes)
        self.moltypes = set(moltypes)

        num_vals = set(num_vals)
        if 0 in num_vals:
            num_vals.remove(0)
        self.num_vals = num_vals

        scaled_vals = set(scaled_vals)
        if 0 in scaled_vals:
            scaled_vals.remove(0)
        self.scaled_vals = scaled_vals

    def __repr__(self):
        return "SignatureParams({}, {}, {}, {})".format(repr(self.ksizes),
                                                        repr(self.moltypes),
                                                        repr(self.num_vals),
                                                        repr(self.scaled_vals))

    def select_ksize(self, ksize):
        if ksize in self.ksizes:
            self.ksizes = { ksize }       # this one only!
            return True
        return False

    def intersect_ksizes(self, ksizes):
        assert ksizes.intersection(self.ksizes)
        self.ksizes.intersection_update(ksizes)

    def select_moltype(self, moltype):
        if moltype in self.moltypes:
            self.moltypes = { moltype }
            return True
        return False

    def intersect_moltypes(self, moltypes):
        assert moltypes.intersection(self.moltypes)
        self.moltypes.intersection_update(moltypes)

    def contains_compatible(self, ksize, moltype, is_scaled):
        if ksize in self.ksizes and moltype in self.moltypes:
            if is_scaled and self.scaled_vals:
                return True
            elif not is_scaled and self.num_vals:
                return True

        return False


def load_SBT_with_params(filename):
    "Load an SBT, construct a params object, return both."
    try:
        tree = SBT.load(filename, leaf_loader=SigLeaf.load)
        num_val, scaled_val, ksize, moltype = get_SBT_info(tree)
    except (ValueError, EnvironmentError):
        return None, None

    assert not (num_val and scaled_val)
    params = SignatureParams([ksize], [moltype], [num_val], [scaled_val])

    return tree, params


def load_LCA_with_params(filename):
    "Load an LCA database, construct a params object, return both."
    try:
        lca_db = lca_utils.LCA_Database()
        lca_db.load(filename)
    except (ValueError, TypeError, EnvironmentError):
        return None, None

    ksize = lca_db.ksize
    # LCA databases are always DNA
    moltype = 'DNA'
    # LCA databases are always made from --scaled signatures
    scaled = lca_db.scaled

    params = SignatureParams([ksize], [moltype], {}, [scaled])

    return lca_db, params


def load_signatures_from_directory_with_params(path):
    """Traverse into a directory path and load all signatures; return w/params.

    This returns union of params across all signatures.
    """
    if not os.path.isdir(path):
        return None, None

    ksizes = set()
    moltypes = set()
    scaled_vals = set()
    num_vals = set()

    siglist = []
    for sigfile in traverse_find_sigs([path]):
        try:
            sigs = sig.load_signatures(sigfile)
        except Exception:        # ignore errors!
            continue

        for ss in sigs:
            # construct union of params across all sigs
            ksizes.add(ss.minhash.ksize)
            moltypes.add(get_moltype(ss))
            scaled_vals.add(ss.minhash.scaled)
            num_vals.add(ss.minhash.num)

            siglist.append(ss)

    if not siglist:
        return None, None

    linear_index = LinearIndex(siglist, path)
    params = SignatureParams(ksizes, moltypes, num_vals, scaled_vals)

    return linear_index, params


def load_signatures_from_file_with_params(filename):
    assert not os.path.isdir(filename)

    ksizes = set()
    moltypes = set()
    scaled_vals = set()
    num_vals = set()

    siglist = []
    if not os.path.exists(filename):
        error("Cannot open file {}", filename)
        raise IOError

    sigs = sig.load_signatures(filename)

    for ss in sigs:
        # construct union of params across all sigs
        ksizes.add(ss.minhash.ksize)
        moltypes.add(get_moltype(ss))
        scaled_vals.add(ss.minhash.scaled)
        num_vals.add(ss.minhash.num)

        siglist.append(ss)

    if not siglist:
        return None, None

    linear_index = LinearIndex(siglist, filename)
    params = SignatureParams(ksizes, moltypes, num_vals, scaled_vals)

    return linear_index, params


def load_target_with_params(filename):
    ordered_loaders = [load_signatures_from_directory_with_params,
                       load_SBT_with_params,
                       load_LCA_with_params,
                       load_signatures_from_file_with_params]

    for loader in ordered_loaders:
        target, params = loader(filename)
        if target:
            return target, params

    return None, None


def filter_compatible_signatures(query, siglist, force=False):
    for ss in siglist:
        if check_signatures_are_compatible(query, ss):
            yield ss
        else:
            if not force:
                raise ValueError("incompatible signature")


# @CTB only for scaled/num values !?
def check_signatures_are_compatible(query, subject):
    if query.minhash.ksize != subject.minhash.ksize:
        return 0

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

def get_SBT_info(tree):
    # get a minhash from the tree
    leaf = next(iter(tree.leaves()))
    tree_ss = leaf.data
    tree_mh = tree_ss.minhash

    moltype = get_moltype(tree_ss)

    return tree_mh.num, tree_mh.scaled, tree_mh.ksize, moltype

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


class SearchDBLoader2(object):
    def __init__(self, require_scaled):
        self.required_scaled = require_scaled
        self.is_query_loaded = False
        self.is_query_checked = False
        self.query_sigs = None
        self.query_params = None
        self.query_filename = None
        self.chosen_query = None

        self.is_args_selector_loaded = False
        self.moltype_selector = None
        self.ksize_selector = None

        self.database_params = []

    def load_query(self, query_filename):
        "Load all the signatures in the given filename, as potential queries"
        assert not self.is_query_loaded
        try:
            (query_sigs, query_params) = load_signatures_from_file_with_params(query_filename)
        except IOError as e:
            query_sigs = None

        if not query_sigs:
            return False

        self.query_sigs = query_sigs
        self.query_params = query_params
        self.is_query_loaded = True
        return True

    def parse_args_selectors(self, args):
        "Parse an ArgumentParser instance for moltype & k-mer size."
        assert not self.is_args_selector_loaded
        moltype = calculate_moltype(args)   # rename this function!!
        ksize = None
        if args.ksize:
            ksize = args.ksize

        self.moltype_selector = moltype
        self.ksize_selector = ksize
        self.is_args_selector_loaded = True

    def check_query_against_arg_selectors(self):
        "Narrow down query parameters against arguments; potentially set 'em."
        if not self.is_args_selector_loaded:
            raise Exception
        if not self.is_query_loaded:
            raise Exception

        assert not self.is_query_checked

        moltype_ok = True
        if self.moltype_selector:
            if not self.query_params.select_moltype(self.moltype_selector):
                moltype_ok = False        # fail! not compatible.
            
        ksize_ok = True
        if self.ksize_selector:
            if not self.query_params.select_ksize(self.ksize_selector):
                ksize_ok = False          # fail! not compatible

        self.is_query_checked = True
        return moltype_ok and ksize_ok

    def add_database(self, identifier, params):
        if not self.is_query_checked:
            raise Exception

        ksize_intersection = params.ksizes.intersection(self.query_params.ksizes)
        moltype_intersection = params.moltypes.intersection(self.query_params.moltypes)

        if len(ksize_intersection) and len(moltype_intersection):
            # save it!
            self.database_params.append((identifier, params))

            # can we nail it down exactly?
            if len(ksize_intersection) == 1 and len(moltype_intersection) == 1:
                self.query_params.select_ksize(ksize_intersection.pop())
                self.query_params.select_moltype(moltype_intersection.pop())
            else:
                # narrow down the query some more.
                self.query_params.intersect_moltypes(params.moltypes)
                self.query_params.intersect_ksizes(params.ksizes)

            return True

        return False

    def decide_query(self):
        if len(self.query_params.ksizes) == 1 and len(self.query_params.moltypes) == 1:
            ksize = next(iter(self.query_params.ksizes))
            moltype = next(iter(self.query_params.moltypes))

            siglist = []
            for sig in self.query_sigs.signatures():
                if sig.minhash.ksize == ksize and get_moltype(sig) == moltype:
                    siglist.append(sig)

            if len(siglist) == 1:
                self.chosen_query = siglist[0]
                return True
#        else:
#            print('XXX', len(self.query_params.ksizes), len(self.query_params.moltypes))

        return False
        

class SearchDatabaseLoader(object):
    def __init__(self, filenames, is_similarity_query, traverse):
        self.filenames = filenames
        self.is_similarity_query = is_similarity_query
        self.traverse = traverse
        self.databases = []

        self.ksizes = set()
        self.moltypes = set()
        self.num_db = set()
        self.scaled_db = set()

    @property
    def ksize(self):
        if not self.ksizes:
            raise ValueError('no ksizes yet')
        if len(self.ksizes) > 1:
            raise ValueError('too many ksizes!')
        return next(iter(self.ksizes))

    @ksize.setter
    def ksize(self, value):
        self.ksizes.add(value)
        if len(self.ksizes) > 1:
            raise ValueError('too many ksizes!')

    @property
    def moltype(self):
        if not self.moltypes:
            raise ValueError('no moltypes yet')
        if len(self.moltypes) > 1:
            raise ValueError('too many moltypes!')
        return next(iter(self.moltypes))

    @moltype.setter
    def moltype(self, value):
        self.moltypes.add(value)
        if len(self.moltypes) > 1:
            raise ValueError('too many moltypes')

    @property
    def scaled(self):
        if self.scaled_db and self.num_db:
            raise ValueError('cannot have both scaled and num DBs')
        if not self.scaled_db and not self.num_db:
            raise ValueError('no information yet on scaled/num')
        if self.scaled_db:
            assert not self.num_db
            return True
        else: # num_db
            return False

    def load_list_of_signatures(self, filename):
        "Load a directory or path full of signatures."
        # are we collecting signatures from a directory/path?
        if not (self.traverse and os.path.isdir(filename)):
            return False

        for sigfile in traverse_find_sigs([filename]):
            try:
                siglist = sig.load_signatures(sigfile)
                linear = LinearIndex(siglist, filename=sigfile)
                self.databases.append((linear, filename, 'linear'))
                notify('loaded {} signatures from {}', len(linear),
                       sigfile, end='\r')
                n_signatures += len(linear)
            except Exception:        # ignore errors with traverse
                pass

        return True

    def load_SBT(self, filename):
        "Load an SBT database."
        try:
            tree = SBT.load(filename, leaf_loader=SigLeaf.load)

            is_num, is_scaled, ksize, moltype = get_SBT_info(tree)
            if is_num:
                self.num_db.add(filename)
            if is_scaled:
                self.scaled_db.add(filename)
            self.ksizes.add(ksize)
            self.moltypes.add(moltype)

            self.databases.append((tree, filename, 'SBT'))
            notify('loaded SBT {}', filename, end='\r')
        except (ValueError, EnvironmentError):
            # not an SBT
            return False

        return True

    def load_LCA(self, filename):
        "Load an LCA database."
        try:
            lca_db = lca_utils.LCA_Database()
            lca_db.load(filename)

            self.ksizes.add(lca_db.ksize)
            # LCA databases are always made from --scaled signatures
            self.scaled_db.add(filename)
            # LCA databases are always DNA
            self.moltypes.add('DNA')

            notify('loaded LCA {}', filename, end='\r')

            self.databases.append((lca_db, filename, 'LCA'))
        except (ValueError, TypeError, EnvironmentError):
            # not an LCA database
            return False

        return True

    def load_signature_file(self, filename):
        "Load a one or more signatures from a file."
        try:
            siglist = sig.load_signatures(filename)
            siglist = list(siglist)
            if len(siglist) == 0:         # file not found, or parse error?
                raise ValueError

            linear = LinearIndex(siglist, filename=filename)
            self.databases.append((linear, filename, 'linear'))

            notify('loaded {} signatures from {}', len(linear), filename,
                   end='\r')
        except (EnvironmentError, ValueError):
            return False

        return True

    def load_all(self):
        for filename in self.filenames:
            notify('loading from {}...', filename, end='\r')

            loaded = self.load_list_of_signatures(filename)
            if not loaded:
                loaded = self.load_SBT(filename)
            if not loaded:
                loaded = self.load_LCA(filename)
            if not loaded:
                loaded = self.load_signature_file(filename)

            if not loaded:
                error("\nCannot load signatures from file '{}'", filename)
                sys.exit(-1)

            # if we loaded a database, and it is incompatible, break out here!
            if self.num_db and self.scaled_db:
                if filename in self.num_db:
                    assert self.scaled_db
                    assert len(self.num_db) == 1
                    error('file {} is incompatible with the other databases',
                          filename)
                    error('the signatures in it were calculated with --num')
                    error('the other databases have --scaled signatures')
                else:
                    assert filename in self.scaled_db
                    assert self.num_db
                    assert len(self.scaled_db) == 1
                    error('file {} is incompatible with the other databases',
                          filename)
                    error('the signatures in it were calculated with --scaled')
                    error('the other databases have --num signatures')

    def filter_signatures(self):
        """
        Run through all the signatures that we loaded, filtering out the
        ones that are incompatible with set parameters.
        """
        ksize = self.ksize
        moltype = self.moltype
        is_scaled = self.scaled

        def select_matching_sigs(sigiter):
            for sig in sigiter:
                # eliminate for a variety of reasons...
                if sig.minhash.ksize != ksize:
                    continue
                if get_moltype(sig) != moltype:
                    continue
                if sig.minhash.scaled and not is_scaled:
                    continue
                if not sig.minhash.scaled and is_scaled:
                    continue

                # keep!
                yield sig

        for (obj, filename, objtype) in self.databases:
            if objtype == 'linear':
                siglist = list(select_matching_sigs(obj.signatures()))

                if not siglist:
                    error("No compatible signatures in {}", filename)
                    sys.exit(-1)

                # @CTB we shouldn't reach deeply into the guts here :)
                obj._signatures = siglist

    def summarize_files(self):
        "Summarize counts of signatures and databases"
        n_signatures = 0
        n_databases = 0

        for (obj, filename, objtype) in self.databases:
            if objtype == 'linear':
                n_signatures += len(obj)
            else:
                assert objtype in ('SBT', 'LCA')
                n_databases += 1

        return n_signatures, n_databases


def load_dbs_and_sigs(filenames, query, is_similarity_query, traverse=False):
    """
    Load one or more SBTs, LCAs, and/or signatures.

    Check for compatibility with query.
    """
    loader = SearchDatabaseLoader(filenames, is_similarity_query, traverse)

    # this loads all the databases and signatures without filtering
    loader.load_all()

    # now we need to figure out what the parameters are ...
    ksize = query.minhash.ksize
    scaled = query.minhash.scaled
    moltype = get_moltype(query)

    # did we load any databases? if so, check against query values. @CTB
    if loader.ksizes:
        assert loader.moltypes
        assert loader.scaled_db or loader.num_db

        if ksize not in loader.ksizes:
            assert 0
        if moltype not in loader.moltypes:
            assert 0
        if scaled and not loader.scaled_db:
            assert 0
        if not scaled and not loader.num_db:
            assert 0
    else:                                 # no databases loaded. @CTB
        loader.ksizes.add(ksize)
        loader.moltypes.add(moltype)
        if scaled:
            loader.scaled_db.add('foo')
        else:
            loader.num_db.add('foo')

    # ok, see if we can filter the signatures down to match.
    loader.filter_signatures()

    # get final counts --
    n_signatures, n_databases = loader.summarize_files()

    notify(' '*79, end='\r')
    if n_signatures and n_databases:
        notify('loaded {} signatures and {} databases total.',
                    n_signatures, n_databases)
    elif n_signatures:
        notify('loaded {} signatures.', n_signatures)
    elif n_databases:
        notify('loaded {} databases.', n_databases)
    else:
        sys.exit(-1)

    return loader.databases


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
