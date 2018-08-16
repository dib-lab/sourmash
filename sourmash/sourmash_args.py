"""
Utility functions for dealing with input args to the sourmash command line.
"""
import collections
import sys
import os
from . import signature
from .logging import notify, error

from . import signature as sig
from .sbt import SBT, Node
from .sbtmh import SigLeaf, select_signature

DEFAULT_LOAD_K=31


def add_moltype_args(parser):
    parser.add_argument('--protein', dest='protein', action='store_true',
                        help='choose a protein signature (default: False)')
    parser.add_argument('--no-protein', dest='protein',
                        action='store_false',
                        help='do not choose a protein signature')
    parser.set_defaults(protein=False)

    parser.add_argument('--dna', dest='dna', default=None,
                        action='store_true',
                        help='choose a DNA signature (default: True)')
    parser.add_argument('--no-dna', dest='dna', action='store_false',
                        help='do not choose a DNA signature')
    parser.set_defaults(dna=None)


def add_construct_moltype_args(parser):
    parser.add_argument('--protein', dest='protein', action='store_true',
                        help='build protein signatures (default: False)')
    parser.add_argument('--no-protein', dest='protein',
                        action='store_false',
                        help='do not build protein signatures')
    parser.set_defaults(protein=False)

    parser.add_argument('--dna', dest='dna', default=None,
                        action='store_true',
                        help='build DNA signatures (default: True)')
    parser.add_argument('--no-dna', dest='dna', action='store_false',
                        help='do not build DNA signatures')
    parser.set_defaults(dna=True)


def add_ksize_arg(parser, default):
    parser.add_argument('-k', '--ksize', default=None, type=int,
                        help='k-mer size (default: {d})'.format(d=default))


def get_moltype(sig, require=False):
    if sig.minhash.is_molecule_type('DNA'):
        moltype = 'DNA'
    elif sig.minhash.is_molecule_type('protein'):
        moltype = 'protein'
    else:
        raise ValueError('unknown molecule type for sig {}'.format(sig.name()))

    return moltype


def calculate_moltype(args, default=None):
    if args.protein:
        if args.dna is True:
            error('cannot specify both --dna and --protein!')
            sys.exit(-1)
        args.dna = False

    moltype = default
    if args.protein:
        moltype = 'protein'
    elif args.dna:
        moltype = 'DNA'

    return moltype


def load_query_signature(filename, ksize, select_moltype):
    try:
        sl = signature.load_signatures(filename,
                                       ksize=ksize,
                                       select_moltype=select_moltype,
                                       do_raise=True)
        sl = list(sl)
    except IOError:
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
        error('need exactly one. Specify --ksize or --dna/--protein.')
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

                yield filename, query, query_moltype, query_ksize

            if len(self.ksizes) > 1 or len(self.moltypes) > 1:
                raise ValueError('multiple k-mer sizes/molecule types present')


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
    query_mh = query.minhash

    leaf = next(iter(tree.leaves()))
    tree_mh = select_signature(leaf, query)
    if tree_mh is None:
        error("tree '{}' is not compatible", treename)
        error('this is different from query ksize of {}.', query_mh.ksize)
        return 0

    tree_mh = tree_mh.minhash

    tree_ksize = get_ksize(tree)
    if tree_ksize != query_mh.ksize:
        error("ksize on tree '{}' is {};", treename, tree_ksize)
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


def get_ksize(tree):
    """Walk nodes in `tree` to find out ksizes"""
    try:
        return tree.nodes[0].data.ksize()
    except AttributeError:
        return get_all_ksizes(tree)


def get_all_ksizes(tree):
    """Walk nodes in `tree` to find out all ksizes"""
    ksizes = set()
    for node in tree.nodes.values():
        if isinstance(node, Node):
            return node.data.ksize()
        elif isinstance(node, SigLeaf):
            if isinstance(node.data, collections.Sequence):
                ksizes |= {q.minhash.ksize for q in node.data}
            else:
                ksizes |= node.data.minhash.ksize

    if len(ksizes) == 1:
        return ksizes.pop()

    return ksizes


def load_sbts_and_sigs(filenames, query, is_similarity_query, traverse=False):
    query_ksize = query.minhash.ksize
    query_moltype = get_moltype(query)

    n_signatures = 0
    n_databases = 0
    databases = []
    for sbt_or_sigfile in filenames:
        if traverse and os.path.isdir(sbt_or_sigfile):
            for sigfile in traverse_find_sigs([sbt_or_sigfile]):
                try:
                    siglist = sig.load_signatures(sigfile,
                                                  ksize=query_ksize,
                                                  select_moltype=query_moltype)
                    siglist = filter_compatible_signatures(query, siglist, 1)
                    siglist = list(siglist)
                    databases.append((siglist, sbt_or_sigfile, False))
                    notify('loaded {} signatures from {}', len(siglist),
                           sigfile, end='\r')
                    n_signatures += len(siglist)
                except:                       # ignore errors with traverse
                    pass

            # done! jump to beginning of main 'for' loop
            continue

        # no traverse? try loading as an SBT.
        try:
            tree = SBT.load(sbt_or_sigfile, leaf_loader=SigLeaf.load)

            if not check_tree_is_compatible(sbt_or_sigfile, tree, query,
                                            is_similarity_query):
                sys.exit(-1)

            databases.append((tree, sbt_or_sigfile, True))
            notify('loaded SBT {}', sbt_or_sigfile, end='\r')
            n_databases += 1

            # done! jump to beginning of main 'for' loop
            continue
        except (ValueError, EnvironmentError):
            # not an SBT - try as a .sig
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
            siglist = list(siglist)

            databases.append((siglist, sbt_or_sigfile, False))
            notify('loaded {} signatures from {}', len(siglist),
                   sbt_or_sigfile, end='\r')
            n_signatures += len(siglist)
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
