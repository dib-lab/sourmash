"""
Utility functions for lowest-common-ancestor analysis tools.
"""
from __future__ import print_function, division
import sys
import json
import gzip
from os.path import exists
from collections import OrderedDict, namedtuple, defaultdict, Counter
import functools


__all__ = ['taxlist', 'zip_lineage', 'build_tree', 'find_lca',
           'load_single_database', 'load_databases', 'gather_assignments',
           'count_lca_for_assignments', 'LineagePair', 'display_lineage']

try:                                      # py2/py3 compat
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest

from .._minhash import get_max_hash_for_scaled
from ..logging import notify, error, debug
from ..index import Index

# type to store an element in a taxonomic lineage
LineagePair = namedtuple('LineagePair', ['rank', 'name'])


def cached_property(fun):
    """A memoize decorator for class properties."""
    @functools.wraps(fun)
    def get(self):
        try:
            return self._cache[fun]
        except AttributeError:
            self._cache = {}
        except KeyError:
            pass
        ret = self._cache[fun] = fun(self)
        return ret
    return property(get)


def check_files_exist(*files):
    ret = True
    not_found = []
    for f in files:
        if not exists(f):
            not_found.append(f)
            ret = False

    if len(not_found):
        error('Error! Could not find the following files.'
              ' Make sure the file paths are specified correctly.\n{}'.format('\n'.join(not_found)))

    return ret


# ordered list of taxonomic ranks
def taxlist(include_strain=True):
    """
    Provide an ordered list of taxonomic ranks.
    """
    for k in ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus',
              'species']:
        yield k
    if include_strain:
        yield 'strain'


# produce an ordered list of tax names from lineage
def zip_lineage(lineage, include_strain=True, truncate_empty=False):
    """
    Given an iterable of LineagePair objects, return list of lineage names.

    This utility function handles species/strain and empty lineage entries
    gracefully.

    >>> x = [ LineagePair('superkingdom', 'a'), LineagePair('phylum', 'b') ]
    >>> zip_lineage(x)
    ['a', 'b', '', '', '', '', '', '']

    >>> x = [ LineagePair('superkingdom', 'a'), LineagePair(None, ''), LineagePair('class', 'c') ]
    >>> zip_lineage(x)
    ['a', '', 'c', '', '', '', '', '']
    """

    empty = LineagePair(None, '')

    pairs = zip_longest(taxlist(include_strain=include_strain),
                        lineage, fillvalue=empty)
    pairs = list(pairs)

    # eliminate empty if so requested
    if truncate_empty:
        last_lineage_tup = pairs[-1][1]
        while pairs and last_lineage_tup == empty:
            pairs.pop(-1)
            if pairs:
                last_lineage_tup = pairs[-1][1]

    row = []
    for taxrank, lineage_tup in pairs:
        # validate non-empty tax, e.g. superkingdom/phylum/class in order.
        if lineage_tup != empty and lineage_tup.rank != taxrank:
            raise ValueError('incomplete lineage at {} - is {} instead'.format(taxrank, lineage_tup.rank))

        row.append(lineage_tup.name)
    return row


def display_lineage(lineage, include_strain=True, truncate_empty=True):
    return ";".join(zip_lineage(lineage,
                                include_strain=include_strain,
                                truncate_empty=truncate_empty))


# filter function toreplace blank/na/null with 'unassigned'
filter_null = lambda x: 'unassigned' if x.strip() in \
  ('[Blank]', 'na', 'null', '') else x
null_names = set(['[Blank]', 'na', 'null'])


def build_tree(assignments, initial=None):
    """
    Builds a tree of dictionaries from lists of LineagePair objects
    in 'assignments'.  This tree can then be used to find lowest common
    ancestor agreements/confusion.
    """
    if initial is None:
        tree = {}
    else:
        tree = initial

    if not assignments:
        raise ValueError("empty assignment passed to build_tree")

    for assignment in assignments:
        node = tree

        for lineage_tup in assignment:
            if lineage_tup.name:
                child = node.get(lineage_tup, {})
                node[lineage_tup] = child

                # shift -> down in tree
                node = child

    return tree


def find_lca(tree):
    """
    Given a tree produced by 'find_tree', find the first node with multiple
    children, OR the only leaf in the tree.  Return (lineage_tup, reason),
    where 'reason' is the number of children of the returned node, i.e.e
    0 if it's a leaf and > 1 if it's an internal node.
    """

    node = tree
    lineage = []
    while 1:
        if len(node) == 1:                # descend to only child; track path
            lineage_tup = next(iter(node.keys()))
            lineage.append(lineage_tup)
            node = node[lineage_tup]
        elif len(node) == 0:              # at leaf; end
            return tuple(lineage), 0
        else:                             # len(node) > 1 => confusion!!
            return tuple(lineage), len(node)


class LCA_Database(Index):
    """
    Wrapper class for taxonomic database.

    obj.ident_to_idx: key 'identifier' to 'idx'
    obj.idx_to_lid: key 'idx' to 'lid'
    obj.lid_to_lineage: key 'lid' to tuple of LineagePair objects
    obj.hashval_to_idx: key 'hashval' => set('idx')
    """
    def __init__(self):
        self.ksize = None
        self.scaled = None
        
        self.ident_to_idx = None
        self.idx_to_lid = None
        self.lid_to_lineage = None
        self.hashval_to_idx = None
    
        self.filename = None

    def __repr__(self):
        return "LCA_Database('{}')".format(self.filename)

    def signatures(self):
        from .. import SourmashSignature
        for v in self._signatures.values():
            yield SourmashSignature(v)

    def select(self, ksize=None, moltype=None):
        ok = True
        if ksize is not None and self.ksize != ksize:
            ok = False
        if moltype is not None and moltype != 'DNA':
            ok = False

        if ok:
            return self

        raise ValueError("cannot select LCA on ksize {} / moltype {}".format(ksize, moltype))

    def load(self, db_name):
        "Load from a JSON file."
        xopen = open
        if db_name.endswith('.gz'):
            xopen = gzip.open

        with xopen(db_name, 'rt') as fp:
            load_d = {}
            try:
                load_d = json.load(fp)
            except json.decoder.JSONDecodeError:
                pass

            if not load_d:
                raise ValueError("cannot parse database file '{}' as JSON; invalid format.")

            version = None
            db_type = None
            try:
                version = load_d.get('version')
                db_type = load_d.get('type')
            except AttributeError:
                pass

            if db_type != 'sourmash_lca':
                raise ValueError("database file '{}' is not an LCA db.".format(db_name))

            if version != '2.0' or 'lid_to_lineage' not in load_d:
                raise ValueError("Error! This is an old-style LCA DB. You'll need to build or download a newer one.")

            ksize = int(load_d['ksize'])
            scaled = int(load_d['scaled'])
            self.ksize = ksize
            self.scaled = scaled

            # convert lineage_dict to proper lineages (tuples of LineagePairs)
            lid_to_lineage_2 = load_d['lid_to_lineage']
            lid_to_lineage = {}
            for k, v in lid_to_lineage_2.items():
                v = dict(v)
                vv = []
                for rank in taxlist():
                    name = v.get(rank, '')
                    vv.append(LineagePair(rank, name))

                lid_to_lineage[int(k)] = tuple(vv)
            self.lid_to_lineage = lid_to_lineage

            # convert hashval -> lineage index keys to integers (looks like
            # JSON doesn't have a 64 bit type so stores them as strings)
            hashval_to_idx_2 = load_d['hashval_to_idx']
            hashval_to_idx = {}

            for k, v in hashval_to_idx_2.items():
                hashval_to_idx[int(k)] = v
            self.hashval_to_idx = hashval_to_idx

            self.ident_to_name = load_d['ident_to_name']
            self.ident_to_idx = load_d['ident_to_idx']

            self.idx_to_lid = {}
            for k, v in load_d['idx_to_lid'].items():
                self.idx_to_lid[int(k)] = v

        self.filename = db_name

    def save(self, db_name):
        "Save to a JSON file."
        xopen = open
        if db_name.endswith('.gz'):
            xopen = gzip.open

        with xopen(db_name, 'wt') as fp:
            # use an OrderedDict to preserve output order
            save_d = OrderedDict()
            save_d['version'] = '2.0'
            save_d['type'] = 'sourmash_lca'
            save_d['license'] = 'CC0'
            save_d['ksize'] = self.ksize
            save_d['scaled'] = self.scaled

            # convert lineage internals from tuples to dictionaries
            d = OrderedDict()
            for k, v in self.lid_to_lineage.items():
                d[k] = dict([ (vv.rank, vv.name) for vv in v ])
            save_d['lid_to_lineage'] = d

            # convert values from sets to lists, so that JSON knows how to save
            save_d['hashval_to_idx'] = \
               dict((k, list(v)) for (k, v) in self.hashval_to_idx.items())

            save_d['ident_to_name'] = self.ident_to_name
            save_d['ident_to_idx'] = self.ident_to_idx
            save_d['idx_to_lid'] = self.idx_to_lid
            save_d['lid_to_lineage'] = self.lid_to_lineage
            
            json.dump(save_d, fp)

    def search(self, query, *args, **kwargs):
        # check arguments
        if 'threshold' not in kwargs:
            raise TypeError("'search' requires 'threshold'")
        threshold = kwargs['threshold']
        do_containment = kwargs.get('do_containment', False)
        ignore_abundance = kwargs.get('ignore_abundance', False)
        mh = query.minhash
        if ignore_abundance:
            mh.track_abundance = False

        results = []
        for x in self.find_signatures(mh, threshold, do_containment):
            (score, match, filename) = x
            results.append((score, match, filename))

        results.sort(key=lambda x: -x[0])
        return results

    def gather(self, query, *args, **kwargs):
        if not query.minhash:
            return []

        results = []
        threshold_bp = kwargs.get('threshold_bp', 0.0)
        threshold = threshold_bp / (len(query.minhash) * self.scaled)
        for x in self.find_signatures(query.minhash, threshold,
                                      containment=True, ignore_scaled=True):
            (score, match, filename) = x
            if score:
                results.append((score, match, filename))
                break

        return results

    def insert(self, node):
        raise NotImplementedError

    def find(self, search_fn, *args, **kwargs):
        raise NotImplementedError

    def downsample_scaled(self, scaled):
        """
        Downsample to the provided scaled value, i.e. eliminate all hashes
        that don't fall in the required range.

        NOTE: we probably need to invalidate some of the dynamically
        calculated members of this object, like _signatures, when we do this.
        But we aren't going to right now.
        """
        if scaled == self.scaled:
            return
        elif scaled < self.scaled:
            raise ValueError("cannot decrease scaled from {} to {}".format(self.scaled, scaled))

        max_hash = get_max_hash_for_scaled(scaled)
        new_hashvals = {}
        for k, v in self.hashval_to_idx.items():
            if k < max_hash:
                new_hashvals[k] = v
        self.hashval_to_idx = new_hashvals
        self.scaled = scaled

    def get_lineage_assignments(self, hashval):
        """
        Get a list of lineages for this hashval.
        """
        x = []

        idx_list = self.hashval_to_idx.get(hashval, [])
        for idx in idx_list:
            lid = self.idx_to_lid.get(idx, None)
            if lid is not None:
                lineage = self.lid_to_lineage[lid]
                x.append(lineage)

        return x

    @cached_property
    def _signatures(self):
        "Create a _signatures member dictionary that contains {idx: minhash}."
        from .. import MinHash

        minhash = MinHash(n=0, ksize=self.ksize, scaled=self.scaled)

        debug('creating signatures for LCA DB...')
        sigd = defaultdict(minhash.copy_and_clear)
        temp_vals = defaultdict(list)

        for (k, v) in self.hashval_to_idx.items():
            for vv in v:
                temp_hashes = temp_vals[vv]
                temp_hashes.append(k)

                # 50 is an arbitrary number. If you really want
                # to micro-optimize, list is resized and grow in this pattern:
                # 0, 4, 8, 16, 25, 35, 46, 58, 72, 88, ...
                # (from https://github.com/python/cpython/blob/b2b4a51f7463a0392456f7772f33223e57fa4ccc/Objects/listobject.c#L57)
                if len(temp_hashes) > 50:
                    sigd[vv].add_many(temp_hashes)

                    # Sigh, python 2... when it goes away,
                    # we can do `temp_hashes.clear()` instead.
                    del temp_vals[vv]

        # We loop temp_vals again to add any remainder hashes
        # (each list of hashes is smaller than 50 items)
        for sig, vals in temp_vals.items():
            sigd[sig].add_many(vals)

        debug('=> {} signatures!', len(sigd))
        return sigd

    def find_signatures(self, minhash, threshold, containment=False,
                       ignore_scaled=False):
        """
        Do a Jaccard similarity or containment search.
        """
        # make sure we're looking at the same scaled value as database
        if self.scaled > minhash.scaled:
            minhash = minhash.downsample_scaled(self.scaled)
        elif self.scaled < minhash.scaled and not ignore_scaled:
            # note that containment can be calculated w/o matching scaled.
            raise ValueError("lca db scaled is {} vs query {}; must downsample".format(self.scaled, minhash.scaled))

        query_mins = set(minhash.get_mins())

        # collect matching hashes:
        c = Counter()
        for hashval in query_mins:
            idx_list = self.hashval_to_idx.get(hashval, [])
            for idx in idx_list:
                c[idx] += 1

        debug('number of matching signatures for hashes: {}', len(c))

        for idx, count in c.items():
            ident = self.idx_to_ident[idx]
            name = self.ident_to_name[ident]

            match_mh = self._signatures[idx]
            match_size = len(match_mh)

            debug('count: {}; query_mins: {}; match size: {}',
                  count, len(query_mins), match_size)

            if containment:
                score = count / len(query_mins)
            else:
                score = count / (len(query_mins) + match_size - count)

            debug('score: {} (containment? {}), threshold: {}',
                  score, containment, threshold)

            if score >= threshold:
                from .. import SourmashSignature
                match_sig = SourmashSignature(match_mh, name=name)

                yield score, match_sig, self.filename

    @cached_property
    def lineage_to_lids(self):
        d = defaultdict(set)
        for lid, lineage in self.lid_to_lineage.items():
            d[lineage].add(lid)
        return d

    @cached_property
    def lid_to_idx(self):
        d = defaultdict(set)
        for idx, lid in self.idx_to_lid.items():
            d[lid].add(idx)
        return d

    @cached_property
    def idx_to_ident(self):
        d = defaultdict(set)
        for ident, idx in self.ident_to_idx.items():
            assert idx not in d
            d[idx] = ident
        return d


def load_single_database(filename, verbose=False):
    "Load a single LCA database; return (db, ksize, scaled)"
    dblist, ksize, scaled = load_databases([filename], verbose=verbose)
    return dblist[0], ksize, scaled


def load_databases(filenames, scaled=None, verbose=True):
    "Load multiple LCA databases; return (dblist, ksize, scaled)"
    ksize_vals = set()
    scaled_vals = set()
    dblist = []

    # load all the databases
    for db_name in filenames:
        if verbose:
            notify(u'\r\033[K', end=u'', file=sys.stderr)
            notify('... loading database {}'.format(db_name), end='\r',
                  file=sys.stderr)

        lca_db = LCA_Database()
        lca_db.load(db_name)

        ksize_vals.add(lca_db.ksize)
        if len(ksize_vals) > 1:
            raise Exception('multiple ksizes, quitting')

        if scaled and scaled > lca_db.scaled:
            lca_db.downsample_scaled(scaled)
        scaled_vals.add(lca_db.scaled)

        dblist.append(lca_db)

    ksize = ksize_vals.pop()
    scaled = scaled_vals.pop()

    if verbose:
        notify(u'\r\033[K', end=u'')
        notify('loaded {} LCA databases. ksize={}, scaled={}', len(dblist),
               ksize, scaled)

    return dblist, ksize, scaled


def gather_assignments(hashvals, dblist):
    """
    Gather assignments from across all the databases for all the hashvals.
    """
    assignments = defaultdict(set)
    for hashval in hashvals:
        for lca_db in dblist:
            lineages = lca_db.get_lineage_assignments(hashval)
            if lineages:
                assignments[hashval].update(lineages)

    return assignments


def count_lca_for_assignments(assignments):
    """
    For each hashval, count the LCA across its assignments.
    """
    counts = Counter()
    for hashval in assignments:

        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover lowest-common-ancestor.
        lineages = assignments[hashval]
        tree = build_tree(lineages)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = find_lca(tree)
        counts[lca] += 1

    return counts
