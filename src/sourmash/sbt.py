#!/usr/bin/env python
"""
An implementation of sequence bloom trees, Solomon & Kingsford, 2015.
"""


from collections import namedtuple, Counter
from collections.abc import Mapping

from copy import copy
import json
import math
import os
from random import randint, random
import sys
from tempfile import NamedTemporaryFile
from cachetools import Cache

from .exceptions import IndexNotSupported
from .sbt_storage import FSStorage, IPFSStorage, RedisStorage, ZipStorage
from .logging import error, notify, debug
from .index import Index, IndexSearchResult

from .nodegraph import Nodegraph, extract_nodegraph_info, calc_expected_collisions

STORAGES = {
    'FSStorage': FSStorage,
    'IPFSStorage': IPFSStorage,
    'RedisStorage': RedisStorage,
    'ZipStorage': ZipStorage,
}


NodePos = namedtuple("NodePos", ["pos", "node"])


class GraphFactory(object):
    """Build new nodegraphs (Bloom filters) of a specific (fixed) size.

    Parameters
    ----------
    ksize: int
        k-mer size.
    starting_size: int
        size (in bytes) for each nodegraph table.
    n_tables: int
        number of nodegraph tables to be used.
    """

    def __init__(self, ksize, starting_size, n_tables):
        self.ksize = ksize
        self.starting_size = starting_size
        self.n_tables = n_tables

    def __call__(self):
        return Nodegraph(self.ksize, self.starting_size, self.n_tables)

    def init_args(self):
        return (self.ksize, self.starting_size, self.n_tables)


class _NodesCache(Cache):
    """A cache for SBT nodes that calls .unload() when the node is removed from cache.

    This is adapted from the LFU cache in https://github.com/tkem/cachetools,
    but removing the largest node ids first (those near the bottom/leaves of
    the SBT).
    """

    def __init__(self, maxsize, getsizeof=None):
        Cache.__init__(self, maxsize, getsizeof)
        self.__counter = Counter()

    def __getitem__(self, key, cache_getitem=Cache.__getitem__):
        value = cache_getitem(self, key)
        self.__counter[key] -= 1
        return value

    def __setitem__(self, key, value, cache_setitem=Cache.__setitem__):
        cache_setitem(self, key, value)
        self.__counter[key] -= 1

    def __delitem__(self, key, cache_delitem=Cache.__delitem__):
        cache_delitem(self, key)
        del self.__counter[key]

    def popitem(self):
        """Remove and return the `(key, value)` pair least recently used."""
        try:
            # Select least frequently used keys,
            # limit to 50 items to avoid dealing with huge lists
            common = self.__counter.most_common()[:50]

            # common might include different values, so let's use
            # only keys that have the same value as the first one
            # (all those with the same count are least frequently used items)
            count = common[0][1]

            # we want to remove the item closest to the leaves,
            # and since node ids increase as they get farther from the root
            # we just need to select the maximum key/node id
            (key, _) = max(c for c in common if c[1] == count)
        except IndexError:
            msg = '%s is empty' % self.__class__.__name__
            raise KeyError(msg) from None
        else:
            value = self.pop(key)
            value.unload()
            return (key, value)


class SBT(Index):
    """A Sequence Bloom Tree implementation allowing generic internal nodes and leaves.

    The default node and leaf format is a Bloom Filter (like the original implementation),
    but we also provide a MinHash leaf class (in the sourmash.sbtmh.SigLeaf class)

    Parameters
    ----------
    factory: Factory
        Callable for generating new datastores for internal nodes.
    d: int
        Number of children for each internal node. Defaults to 2 (a binary tree)
    storage: Storage, default: None
        A Storage is any place where we can save and load data for the nodes.
        If set to None, will use a FSStorage.
    cache_size: int, default None
        Number of internal nodes to cache in memory.
        If set to None, will not remove any nodes from memory
        (cache grows without bounds).

    Notes
    -----
    We use two dicts to store the tree structure: One for the internal nodes,
    and another for the leaves (datasets).
    """
    is_database = True

    def __init__(self, factory, *, d=2, storage=None, cache_size=None):
        self.factory = factory
        self._nodes = {}
        self._missing_nodes = set()
        self._leaves = {}
        self.d = d
        self.next_node = 0
        self.storage = storage
        if cache_size is None:
            cache_size = sys.maxsize
        self._nodescache = _NodesCache(maxsize=cache_size)
        self._location = None

    @property
    def location(self):
        return self._location

    def signatures(self):
        for k in self.leaves():
            yield k.data

    def select(self, ksize=None, moltype=None, num=0, scaled=0,
               containment=False):
        """Make sure this database matches the requested requirements.

        Will always raise ValueError if a requirement cannot be met.

        The only tricky bit here is around downsampling: if the scaled
        value being requested is higher than the signatures in the
        SBT, we can use the SBT for containment but not for
        similarity. This is because:

        * if we are doing containment searches, the intermediate nodes
          can still be used for calculating containment of signatures
          with higher scaled values. This is because only hashes that match
          in the higher range are used for containment scores.
        * however, for similarity, _all_ hashes are used, and we cannot
          implicitly downsample or necessarily estimate similarity if
          the scaled values differ.
        """
        # pull out a signature from this collection -
        first_sig = next(iter(self.signatures()))
        db_mh = first_sig.minhash

        # check ksize.
        if ksize is not None and db_mh.ksize != ksize:
            raise ValueError(f"search ksize {ksize} is different from database ksize {db_mh.ksize}")

        # check moltype.
        if moltype is not None and db_mh.moltype != moltype:
            raise ValueError(f"search moltype {moltype} is different from database moltype {db_mh.moltype}")

        # containment requires 'scaled'.
        if containment:
            if not scaled:
                raise ValueError("'containment' requires 'scaled' in SBT.select'")
            if not db_mh.scaled:
                raise ValueError("cannot search this SBT for containment; signatures are not calculated with scaled")

        # 'num' and 'scaled' do not mix.
        if num:
            if not db_mh.num:
                raise ValueError(f"this database was created with 'scaled' MinHash sketches, not 'num'")
            if num != db_mh.num:
                raise ValueError(f"num mismatch for SBT: num={num}, {db_mh.num}")

        if scaled:
            if not db_mh.scaled:
                raise ValueError(f"this database was created with 'num' MinHash sketches, not 'scaled'")

            # we can downsample SBTs for containment operations.
            if scaled > db_mh.scaled and not containment:
                raise ValueError(f"search scaled value {scaled} is less than database scaled value of {db_mh.scaled}")

        return self

    def new_node_pos(self, node):
        if not self._nodes:
            self.next_node = 1
            return 0

        if not self._leaves:
            self.next_node = 2
            return 1

        min_leaf = min(self._leaves.keys())

        next_internal_node = None
        if self.next_node <= min_leaf:
            for i in range(min_leaf):
                if all((i not in self._nodes,
                        i not in self._leaves,
                        i not in self._missing_nodes)):
                    next_internal_node = i
                    break

        if next_internal_node is None:
            self.next_node = max(self._leaves.keys()) + 1
        else:
            self.next_node = next_internal_node

        return self.next_node

    def insert(self, signature):
        "Add a new SourmashSignature in to the SBT."
        from .sbtmh import SigLeaf
        
        leaf = SigLeaf(signature.md5sum(), signature)
        self.add_node(leaf)

    def add_node(self, node):
        pos = self.new_node_pos(node)

        if pos == 0:  # empty tree; initialize w/node.
            n = Node(self.factory, name="internal." + str(pos))
            self._nodes[0] = n
            pos = self.new_node_pos(node)

        # Cases:
        # 1) parent is a Leaf (already covered)
        # 2) parent is a Node (with empty position available)
        #    - add Leaf, update parent
        # 3) parent is a Node (no position available)
        #    - this is covered by case 1
        # 4) parent is None
        #    this can happen with d != 2, in this case create the parent node
        p = self.parent(pos)
        if isinstance(p.node, Leaf):
            # Create a new internal node
            # node and parent are children of new internal node
            n = Node(self.factory, name="internal." + str(p.pos))
            self._nodes[p.pos] = n

            c1, c2 = self.children(p.pos)[:2]

            self._leaves[c1.pos] = p.node
            self._leaves[c2.pos] = node 
            del self._leaves[p.pos]

            for child in (p.node, node):
                child.update(n)
        elif isinstance(p.node, Node):
            self._leaves[pos] = node 
            node.update(p.node)
        elif p.node is None:
            n = Node(self.factory, name="internal." + str(p.pos))
            self._nodes[p.pos] = n
            c1 = self.children(p.pos)[0]
            self._leaves[c1.pos] = node 
            node.update(n)

        # update all parents!
        p = self.parent(p.pos)
        while p:
            self._rebuild_node(p.pos)
            node.update(self._nodes[p.pos])
            p = self.parent(p.pos)

    def _find_nodes(self, search_fn, *args, **kwargs):
        "Search the tree using `search_fn`."

        unload_data = kwargs.get("unload_data", False)

        # initialize search queue with top node of tree
        matches = []
        visited, queue = set(), [0]

        # while the queue is not empty, load each node and apply search
        # function.
        while queue:
            node_p = queue.pop(0)

            # repair while searching.
            node_g = self._leaves.get(node_p, None)
            if node_g is None:
                if node_p in self._nodescache:
                    node_g = self._nodescache[node_p]
                else:
                    node_g = self._nodes.get(node_p, None)
                    if node_g is None:
                        if node_p in self._missing_nodes:
                            self._rebuild_node(node_p)
                            node_g = self._nodes[node_p]
                        else:
                            continue
                    self._nodescache[node_p] = node_g

            # if we have not visited this node before,
            if node_p not in visited:
                visited.add(node_p)

                # apply search fn. If return false, truncate search.
                if search_fn(node_g, *args):

                    # leaf node? it's a match!
                    if isinstance(node_g, Leaf):
                        matches.append(node_g)
                    # internal node? descend.
                    elif isinstance(node_g, Node):
                        if kwargs.get('dfs', True):  # defaults search to dfs
                            for c in self.children(node_p):
                                queue.insert(0, c.pos)
                        else: # bfs
                            queue.extend(c.pos for c in self.children(node_p))

                if unload_data:
                    node_g.unload()

        return matches

    def find(self, search_fn, query, **kwargs):
        """
        Do a Jaccard similarity or containment search, yield results.

        Here 'search_fn' should be an instance of 'JaccardSearch'.

        Queries with higher scaled values than the database
        can still be used for containment search, but not for similarity
        search. See SBT.select(...) for details.
        """
        from .sbtmh import SigLeaf

        search_fn.check_is_compatible(query)

        query_mh = query.minhash

        # figure out downsampling using the first leaf in the tree --
        a_leaf = next(iter(self.leaves()))
        tree_scaled = a_leaf.data.minhash.scaled

        # scaled?
        if tree_scaled:
            assert query_mh.scaled

            # pick the larger scaled of the query & node
            scaled = max(query_mh.scaled, tree_scaled)
            if query_mh.scaled < tree_scaled:
                query_mh = query_mh.downsample(scaled=tree_scaled)

            # provide function to downsample leaf_node as well
            if scaled == tree_scaled:
                downsample_node = lambda x: x
            else:
                def downsample_node(node_mh):
                    return node_mh.downsample(scaled=scaled)
        else:
            assert query_mh.num

            # pick the smaller num of the query & node
            min_num = min(query_mh.num, a_leaf.data.minhash.num)

            # downsample query once:
            if query_mh.num > min_num:
                query_mh = query_mh.downsample(num=min_num)

            # provide function to downsample leaf nodes.
            if min_num == a_leaf.data.minhash.num:
                downsample_node = lambda x: x
            else:
                def downsample_node(node_mh):
                    return node_mh.downsample(num=min_num)

        query_size = len(query_mh)

        # store scores here so we don't need to recalculate
        results = {}

        # construct a function to pass into ._find_nodes; this function
        # will be used to prune tree searches based on internal node scores,
        # in addition to finding leaf nodes.
        def node_search(node, *args, **kwargs):
            is_leaf = False

            # leaf node? downsample so we can do signature comparison.
            if isinstance(node, SigLeaf):
                is_leaf = True

                subj_mh = downsample_node(node.data.minhash)
                subj_size = len(subj_mh)
                subj_mh = subj_mh.flatten()

                assert not subj_mh.track_abundance

                shared_size, total_size = query_mh.intersection_and_union_size(subj_mh)
            else:  # Node / Nodegraph by minhash comparison
                # no downsampling needed --
                shared_size = node.data.matches(query_mh)
                subj_size = node.metadata.get('min_n_below', -1)
                if subj_size == -1:
                    raise ValueError("ERROR: no min_n_below on this tree, cannot search.")
                total_size = subj_size # approximate; do not collect

            # calculate score (exact, if leaf; approximate, if not)
            score = search_fn.score_fn(query_size,
                                       shared_size,
                                       subj_size,
                                       total_size)

            if search_fn.passes(score):
                if is_leaf:     # terminal node? keep.
                    if search_fn.collect(score, node.data):
                        results[node.data] = score
                        return True
                else:           # it's a good internal node, keep.
                    return True

            return False

        # & execute!
        for n in self._find_nodes(node_search, **kwargs):
            yield IndexSearchResult(results[n.data], n.data, self.location)

    def _rebuild_node(self, pos=0):
        """Recursively rebuilds an internal node (if it is not present).

        Parameters
        ----------
        pos: int
            node to be rebuild. Any internal node under it will be rebuild too.
            If you want to rebuild all missing internal nodes you can use pos=0
            (the default).
        """

        node = self._nodes.get(pos, None)
        if node is not None:
            # this node was already build, skip
            return

        node = Node(self.factory, name="internal.{}".format(pos))
        self._nodes[pos] = node
        for c in self.children(pos):
            if c.pos in self._missing_nodes or isinstance(c.node, Leaf):
                cnode = c.node
                if cnode is None:
                    self._rebuild_node(c.pos)
                    cnode = self._nodes[c.pos]
                cnode.update(node)

    def parent(self, pos):
        """Return the parent of the node at position ``pos``.

        If it is the root node (position 0), returns None.

        Parameters
        ----------
        pos: int
            Position of the node in the tree.

        Returns
        -------
        NodePos :
            A NodePos namedtuple with the position and content of the parent node.
        """

        if pos == 0:
            return None
        p = int(math.floor((pos - 1) / self.d))
        if p in self._leaves:
            return NodePos(p, self._leaves[p])

        node = self._nodes.get(p, None)
        return NodePos(p, node)

    def children(self, pos):
        """Return all children nodes for node at position ``pos``.

        Parameters
        ----------
        pos: int
            Position of the node in the tree.

        Returns
        -------
        list of NodePos
            A list of NodePos namedtuples with the position and content of all
            children nodes.
        """
        return [self.child(pos, c) for c in range(self.d)]

    def child(self, parent, pos):
        """Return a child node at position ``pos`` under the ``parent`` node.

        Parameters
        ----------
        parent: int
            Parent node position in the tree.
        pos: int
            Position of the child one under the parent. Ranges from
            [0, arity - 1], where arity is the arity of the SBT
            (usually it is 2, a binary tree).

        Returns
        -------
        NodePos
            A NodePos namedtuple with the position and content of the
            child node.
        """
        cd = self.d * parent + pos + 1
        if cd in self._leaves:
            return NodePos(cd, self._leaves[cd])

        node = self._nodes.get(cd, None)
        return NodePos(cd, node)

    def save(self, path, storage=None, sparseness=0.0, structure_only=False):
        """Saves an SBT description locally and node data to a storage.

        Parameters
        ----------
        path : str
            path to where the SBT description should be saved.
        storage : Storage, optional
            Storage to be used for saving node data.
            Defaults to FSStorage (a hidden directory at the same level of path)
        sparseness : float
            How much of the internal nodes should be saved.
            Defaults to 0.0 (save all internal nodes data),
            can go up to 1.0 (don't save any internal nodes data)
        structure_only: boolean
            Write only the index schema and metadata, but not the data.
            Defaults to False (save data too)

        Returns
        -------
        str
            full path to the new SBT description
        """
        info = {}
        info['d'] = self.d
        info['version'] = 6
        info["index_type"] = self.__class__.__name__  # TODO: check

        # choose between ZipStorage and FS (file system/directory) storage.
        if not path.endswith(".sbt.json"):
            kind = "Zip"
            if not path.endswith('.sbt.zip'):
                path += '.sbt.zip'
            storage = ZipStorage(path)
            backend = "FSStorage"
            name = os.path.basename(path[:-8])
            subdir = '.sbt.{}'.format(name)
            storage_args = FSStorage("", subdir).init_args()
            storage.save(subdir + "/", b"")
            index_filename = os.path.abspath(path)
        else:                             # path.endswith('.sbt.json')
            assert path.endswith('.sbt.json')
            kind = "FS"
            name = os.path.basename(path)
            name = name[:-9]
            index_filename = os.path.abspath(path)

            if storage is None:
                # default storage
                location = os.path.dirname(index_filename)
                subdir = '.sbt.{}'.format(name)

                storage = FSStorage(location, subdir)
                index_filename = os.path.join(location, index_filename)

            backend = [k for (k, v) in STORAGES.items() if v == type(storage)][0]
            storage_args = storage.init_args()

        info['storage'] = {
            'backend': backend,
            'args': storage_args
        }
        info['factory'] = {
            'class': GraphFactory.__name__,
            'args': self.factory.init_args()
        }

        nodes = {}
        leaves = {}
        total_nodes = len(self)
        for n, (i, node) in enumerate(self):
            if node is None:
                continue

            if isinstance(node, Node):
                if random() - sparseness <= 0:
                    continue

            data = {
                # TODO: start using md5sum instead?
                'filename': os.path.basename(node.name),
                'name': node.name
            }

            try:
                node.metadata.pop('max_n_below')
            except (AttributeError, KeyError):
                pass

            data['metadata'] = node.metadata

            if structure_only is False:
                # trigger data loading before saving to the new place
                node.data

                node.storage = storage

                if kind == "Zip":
                    node.save(os.path.join(subdir, data['filename']))
                elif kind == "FS":
                    data['filename'] = node.save(data['filename'])

            if isinstance(node, Node):
                nodes[i] = data
            else:
                leaves[i] = data

            if n % 100 == 0:
                notify("{} of {} nodes saved".format(n+1, total_nodes), end='\r')

        notify("Finished saving nodes, now saving SBT index file.")
        info['nodes'] = nodes
        info['signatures'] = leaves

        if kind == "Zip":
            tree_data = json.dumps(info).encode("utf-8")
            save_path = "{}.sbt.json".format(name)
            storage.save(save_path, tree_data)
            storage.close()

        elif kind == "FS":
            with open(index_filename, 'w') as fp:
                json.dump(info, fp)

        notify("Finished saving SBT index, available at {0}\n".format(index_filename))

        return path


    @classmethod
    def load(cls, location, *, leaf_loader=None, storage=None, print_version_warning=True, cache_size=None):
        """Load an SBT description from a file.

        Parameters
        ----------
        location : str
            path to the SBT description.
        leaf_loader : function, optional
            function to load leaf nodes. Defaults to ``Leaf.load``.
        storage : Storage, optional
            Storage to be used for saving node data.
            Defaults to FSStorage (a hidden directory at the same level of path)

        Returns
        -------
        SBT
            the SBT tree built from the description.
        """
        tempfile = None
        sbt_name = None
        tree_data = None

        if storage is None:
            if ZipStorage.can_open(location):
                storage = ZipStorage(location)
            else:
                if not location.endswith('.sbt.zip'):
                    location2 = location + '.sbt.zip'
                    if ZipStorage.can_open(location2):
                        storage = ZipStorage(location2)

            if storage:
                sbts = storage.list_sbts()
                if len(sbts) == 1:
                    tree_data = storage.load(sbts[0])

                tempfile = NamedTemporaryFile()

                tempfile.write(tree_data)
                tempfile.flush()

                dirname = os.path.dirname(tempfile.name)
                sbt_name = os.path.basename(tempfile.name)

        if sbt_name is None:
            dirname = os.path.dirname(os.path.abspath(location))
            sbt_name = os.path.basename(location)
            if sbt_name.endswith('.sbt.json'):
                sbt_name = sbt_name[:-9]

        sbt_fn = os.path.join(dirname, sbt_name)
        if not sbt_fn.endswith('.sbt.json') and tempfile is None:
            sbt_fn += '.sbt.json'
        with open(sbt_fn) as fp:
            jnodes = json.load(fp)

        if tempfile is not None:
            tempfile.close()

        version = 1
        if isinstance(jnodes, Mapping):
            version = jnodes['version']

        if leaf_loader is None:
            leaf_loader = Leaf.load

        loaders = {
            1: cls._load_v1,
            2: cls._load_v2,
            3: cls._load_v3,
            4: cls._load_v4,
            5: cls._load_v5,
            6: cls._load_v6,
        }

        try:
            loader = loaders[version]
        except KeyError:
            raise IndexNotSupported()

        #if version >= 6:
        #    if jnodes.get("index_type", "SBT") == "LocalizedSBT":
        #        loaders[6] = LocalizedSBT._load_v6

        if version < 3 and storage is None:
            storage = FSStorage(dirname, '.sbt.{}'.format(sbt_name))
        elif storage is None:
            klass = STORAGES[jnodes['storage']['backend']]
            if jnodes['storage']['backend'] == "FSStorage":
                storage = FSStorage(dirname, jnodes['storage']['args']['path'])
            elif storage is None:
                storage = klass(**jnodes['storage']['args'])

        obj = loader(jnodes, leaf_loader, dirname, storage, print_version_warning=print_version_warning, cache_size=cache_size)
        obj._location = location
        return obj

    @staticmethod
    def _load_v1(jnodes, leaf_loader, dirname, storage, *, print_version_warning=True, cache_size=None):

        if jnodes[0] is None:
            raise ValueError("Empty tree!")

        sbt_nodes = {}
        sbt_leaves = {}

        sample_bf = os.path.join(dirname, jnodes[0]['filename'])
        ksize, tablesize, ntables = extract_nodegraph_info(sample_bf)[:3]
        factory = GraphFactory(ksize, tablesize, ntables)

        for i, jnode in enumerate(jnodes):
            if jnode is None:
                continue

            jnode['filename'] = os.path.join(dirname, jnode['filename'])

            if 'internal' in jnode['name']:
                jnode['factory'] = factory
                sbt_node = Node.load(jnode, storage)
                sbt_nodes[i] = sbt_node
            else:
                sbt_node = leaf_loader(jnode, storage)
                sbt_leaves[i] = sbt_node

        tree = SBT(factory, cache_size=cache_size)
        tree._nodes = sbt_nodes
        tree._leaves = sbt_leaves

        return tree

    @classmethod
    def _load_v2(cls, info, leaf_loader, dirname, storage, *, print_version_warning=True, cache_size=None):
        nodes = {int(k): v for (k, v) in info['nodes'].items()}

        if nodes[0] is None:
            raise ValueError("Empty tree!")

        sbt_nodes = {}
        sbt_leaves = {}

        sample_bf = os.path.join(dirname, nodes[0]['filename'])
        k, size, ntables = extract_nodegraph_info(sample_bf)[:3]
        factory = GraphFactory(k, size, ntables)

        for k, node in nodes.items():
            if node is None:
                continue

            node['filename'] = os.path.join(dirname, node['filename'])

            if 'internal' in node['name']:
                node['factory'] = factory
                sbt_node = Node.load(node, storage)
                sbt_nodes[k] = sbt_node
            else:
                sbt_node = leaf_loader(node, storage)
                sbt_leaves[k] = sbt_node

        tree = cls(factory, d=info['d'], cache_size=cache_size)
        tree._nodes = sbt_nodes
        tree._leaves = sbt_leaves

        return tree

    @classmethod
    def _load_v3(cls, info, leaf_loader, dirname, storage, *, print_version_warning=True, cache_size=None):
        nodes = {int(k): v for (k, v) in info['nodes'].items()}

        if not nodes:
            raise ValueError("Empty tree!")

        sbt_nodes = {}
        sbt_leaves = {}

        factory = GraphFactory(*info['factory']['args'])

        max_node = 0
        for k, node in nodes.items():
            if node is None:
                continue

            if 'internal' in node['name']:
                node['factory'] = factory
                sbt_node = Node.load(node, storage)
                sbt_nodes[k] = sbt_node
            else:
                sbt_node = leaf_loader(node, storage)
                sbt_leaves[k] = sbt_node

            max_node = max(max_node, k)

        tree = cls(factory, d=info['d'], storage=storage, cache_size=cache_size)
        tree._nodes = sbt_nodes
        tree._leaves = sbt_leaves
        tree._missing_nodes = {i for i in range(max_node)
                              if i not in sbt_nodes and i not in sbt_leaves}

        if print_version_warning:
            error("WARNING: this is an old index version, please run `sourmash migrate` to update it.")
            error("WARNING: proceeding with execution, but it will take longer to finish!")

        tree._fill_min_n_below()

        return tree

    @classmethod
    def _load_v4(cls, info, leaf_loader, dirname, storage, *, print_version_warning=True, cache_size=None):
        nodes = {int(k): v for (k, v) in info['nodes'].items()}

        if not nodes:
            raise ValueError("Empty tree!")

        sbt_nodes = {}
        sbt_leaves = {}

        factory = GraphFactory(*info['factory']['args'])

        max_node = 0
        for k, node in nodes.items():
            if 'internal' in node['name']:
                node['factory'] = factory
                sbt_node = Node.load(node, storage)
                sbt_nodes[k] = sbt_node
            else:
                sbt_node = leaf_loader(node, storage)
                sbt_leaves[k] = sbt_node

            max_node = max(max_node, k)

        tree = cls(factory, d=info['d'], storage=storage, cache_size=cache_size)
        tree._nodes = sbt_nodes
        tree._leaves = sbt_leaves
        tree._missing_nodes = {i for i in range(max_node)
                              if i not in sbt_nodes and i not in sbt_leaves}

        tree.next_node = max_node

        return tree

    @classmethod
    def _load_v5(cls, info, leaf_loader, dirname, storage, *, print_version_warning=True, cache_size=None):
        nodes = {int(k): v for (k, v) in info['nodes'].items()}
        leaves = {int(k): v for (k, v) in info['leaves'].items()}

        if not leaves:
            raise ValueError("Empty tree!")

        sbt_nodes = {}
        sbt_leaves = {}

        if storage is None:
            klass = STORAGES[info['storage']['backend']]
            if info['storage']['backend'] == "FSStorage":
                storage = FSStorage(dirname, info['storage']['args']['path'])
            elif storage is None:
                storage = klass(**info['storage']['args'])

        factory = GraphFactory(*info['factory']['args'])

        max_node = 0
        for k, node in nodes.items():
            node['factory'] = factory
            sbt_node = Node.load(node, storage)

            sbt_nodes[k] = sbt_node
            max_node = max(max_node, k)

        for k, node in leaves.items():
            sbt_leaf = leaf_loader(node, storage)
            sbt_leaves[k] = sbt_leaf
            max_node = max(max_node, k)

        tree = cls(factory, d=info['d'], storage=storage, cache_size=cache_size)
        tree._nodes = sbt_nodes
        tree._leaves = sbt_leaves
        tree._missing_nodes = {i for i in range(max_node)
                              if i not in sbt_nodes and i not in sbt_leaves}

        return tree

    @classmethod
    def _load_v6(cls, info, leaf_loader, dirname, storage, *, print_version_warning=True, cache_size=None):
        nodes = {int(k): v for (k, v) in info['nodes'].items()}
        leaves = {int(k): v for (k, v) in info['signatures'].items()}

        if not leaves:
            raise ValueError("Empty tree!")

        sbt_nodes = {}
        sbt_leaves = {}

        if storage is None:
            klass = STORAGES[info['storage']['backend']]
            if info['storage']['backend'] == "FSStorage":
                storage = FSStorage(dirname, info['storage']['args']['path'])
            elif storage is None:
                storage = klass(**info['storage']['args'])

        factory = GraphFactory(*info['factory']['args'])

        max_node = 0
        for k, node in nodes.items():
            node['factory'] = factory
            sbt_node = Node.load(node, storage)

            sbt_nodes[k] = sbt_node
            max_node = max(max_node, k)

        for k, node in leaves.items():
            sbt_leaf = leaf_loader(node, storage)
            sbt_leaves[k] = sbt_leaf
            max_node = max(max_node, k)

        tree = cls(factory, d=info['d'], storage=storage, cache_size=cache_size)
        tree._nodes = sbt_nodes
        tree._leaves = sbt_leaves
        tree._missing_nodes = {i for i in range(max_node)
                              if i not in sbt_nodes and i not in sbt_leaves}

        return tree

    def _fill_min_n_below(self):
        """\
        Propagate the smallest hash size below each node up the tree from
        the leaves.
        """
        def fill_min_n_below(node, *args, **kwargs):
            original_min_n_below = node.metadata.get('min_n_below', sys.maxsize)
            min_n_below = original_min_n_below

            children = kwargs['children']
            for child in children:
                if child.node is not None:
                    if isinstance(child.node, Leaf):
                        min_n_below = min(len(child.node.data.minhash), min_n_below)
                    else:
                        child_n = child.node.metadata.get('min_n_below', sys.maxsize)
                        min_n_below = min(child_n, min_n_below)

            if min_n_below == 0:
                min_n_below = 1

            node.metadata['min_n_below'] = min_n_below
            return original_min_n_below != min_n_below

        self._fill_up(fill_min_n_below)

    def _fill_internal(self):

        def fill_nodegraphs(node, *args, **kwargs):
            children = kwargs['children']
            for child in children:
                if child.node is not None:
                    child.node.update(node)
            return True

        self._fill_up(fill_nodegraphs)

    def _fill_up(self, search_fn, *args, **kwargs):
        visited, queue = set(), list(reversed(sorted(self._leaves.keys())))
        debug("started filling up")
        processed = 0
        while queue:
            node_p = queue.pop(0)

            parent = self.parent(node_p)
            if parent is None:
                # we are in the root, no more nodes available to search
                assert len(queue) == 0
                return

            was_missing = False
            if parent.node is None:
                if parent.pos in self._missing_nodes:
                    self._rebuild_node(parent.pos)
                    parent = self.parent(node_p)
                    was_missing = True
                else:
                    continue

            siblings = self.children(parent.pos)

            if node_p not in visited:
                visited.add(node_p)
                for sibling in siblings:
                    visited.add(sibling.pos)
                    try:
                        queue.remove(sibling.pos)
                    except ValueError:
                        pass

                if search_fn(parent.node, children=siblings, *args) or was_missing:
                    queue.append(parent.pos)

            processed += 1
            if processed % 100 == 0:
                debug("processed {}, in queue {}", processed, len(queue), sep='\r')

    def __len__(self):
        internal_nodes = set(self._nodes).union(self._missing_nodes)
        return len(internal_nodes) + len(self._leaves)

    def print_dot(self):
        print("""
        digraph G {
        nodesep=0.3;
        ranksep=0.2;
        margin=0.1;
        node [shape=ellipse];
        edge [arrowsize=0.8];
        """)

        for i, node in self._nodes.items():
            if isinstance(node, Node):
                print('"{}" [shape=box fillcolor=gray style=filled]'.format(
                      node.name))
                for j, child in self.children(i):
                    if child is not None:
                        print('"{}" -> "{}"'.format(node.name, child.name))
        print("}")

    def print(self):
        visited, stack = set(), [0]
        while stack:
            node_p = stack.pop()
            node_g = self._nodes.get(node_p, None)
            if node_p not in visited and node_g is not None:
                visited.add(node_p)
                depth = int(math.floor(math.log(node_p + 1, self.d)))
                print(" " * 4 * depth, node_g)
                if isinstance(node_g, Node):
                    stack.extend(c.pos for c in self.children(node_p)
                                       if c.pos not in visited)

    def __iter__(self):
        for i, node in self._nodes.items():
            yield (i, node)
        for i, node in self._leaves.items():
            yield (i, node)

    def _parents(self, pos=0):
        if pos == 0:
            yield None
        else:
            p = self.parent(pos)
            while p is not None:
                yield p.pos
                p = self.parent(p.pos)

    def leaves(self, with_pos=False):
        for pos, data in self._leaves.items():
            if with_pos:
                yield (pos, data)
            else:
                yield data

    def combine(self, other):
        larger, smaller = self, other
        if len(other) > len(self):
            larger, smaller = other, self

        n = Node(self.factory, name="internal.0", storage=self.storage)
        larger._nodes[0].update(n)
        smaller._nodes[0].update(n)
        new_nodes = {}
        new_nodes[0] = n

        new_leaves = {}

        levels = int(math.ceil(math.log(len(larger), self.d))) + 1
        current_pos = 1
        n_previous = 0
        n_next = 1
        for level in range(1, levels + 1):
            for tree in (larger, smaller):
                for pos in range(n_previous, n_next):
                    if tree._nodes.get(pos, None) is not None:
                        new_node = copy(tree._nodes[pos])
                        new_node.name = "internal.{}".format(current_pos)
                        new_nodes[current_pos] = new_node
                    elif tree._leaves.get(pos, None) is not None:
                        new_node = copy(tree._leaves[pos])
                        new_leaves[current_pos] = new_node
                    current_pos += 1
            n_previous = n_next
            n_next = n_previous + int(self.d ** level)
            current_pos = n_next

        # TODO: do we want to return a new tree, or merge into this one?
        self._nodes = new_nodes
        self._leaves = new_leaves
        return self


class Node(object):
    "Internal node of SBT."

    def __init__(self, factory, name=None, path=None, storage=None):
        self.name = name
        self.storage = storage
        self._factory = factory
        self._data = None
        self._path = path
        self.metadata = dict()

    def __str__(self):
        return '*Node:{name} [occupied: {nb}, fpr: {fpr:.2}]'.format(
                name=self.name, nb=self.data.n_occupied(),
                fpr=calc_expected_collisions(self.data, True, 1.1))

    def save(self, path):
        buf = self.data.to_bytes(compression=1)
        return self.storage.save(path, buf)

    @property
    def data(self):
        if self._data is None:
            if self._path is None:
                self._data = self._factory()
            else:
                data = self.storage.load(self._path)
                self._data = Nodegraph.from_buffer(data)
        return self._data

    @data.setter
    def data(self, new_data):
        self._data = new_data

    def unload(self):
        if self.storage:
            # Don't unload data if there is no Storage
            # TODO: Check that data is actually in the storage?
            self._data = None

    @staticmethod
    def load(info, storage=None):
        new_node = Node(info['factory'],
                        name=info['name'],
                        path=info['filename'],
                        storage=storage)
        new_node.metadata = info.get('metadata', {})
        return new_node

    def update(self, parent):
        parent.data.update(self.data)
        if 'min_n_below' in self.metadata:
            min_n_below = min(parent.metadata.get('min_n_below', sys.maxsize),
                              self.metadata.get('min_n_below'))
            if min_n_below == 0:
                min_n_below = 1
            parent.metadata['min_n_below'] = min_n_below


class Leaf(object):
    def __init__(self, metadata, data=None, name=None, storage=None, path=None):
        self.metadata = metadata

        if name is None:
            name = metadata
        self.name = name

        self.storage = storage

        self._data = data
        self._path = path

    def __str__(self):
        return '**Leaf:{name} [occupied: {nb}, fpr: {fpr:.2}] -> {metadata}'.format(
                name=self.name, metadata=self.metadata,
                nb=self.data.n_occupied(),
                fpr=calc_expected_collisions(self.data, True, 1.1))

    @property
    def data(self):
        if self._data is None:
            data = self.storage.load(self._path)
            self._data = Nodegraph.from_buffer(data)
        return self._data

    @data.setter
    def data(self, new_data):
        self._data = new_data

    def unload(self):
        if self.storage:
            # Don't unload data if there is no Storage
            # TODO: Check that data is actually in the storage?
            self._data = None

    def save(self, path):
        buf = self.data.to_bytes(compression=1)
        return self.storage.save(path, buf)

    def update(self, parent):
        parent.data.update(self.data)

    @classmethod
    def load(cls, info, storage=None):
        return cls(info['metadata'],
                   name=info['name'],
                   path=info['filename'],
                   storage=storage)


def filter_distance(filter_a, filter_b, n=1000):
    """
    Compute a heuristic distance per bit between two Bloom filters.

    Parameters
    ----------
    filter_a : Nodegraph
    filter_b : Nodegraph
    n        : int
        Number of positions to compare (in groups of 8)

    Returns
    -------
    float
        The distance between both filters (from 0.0 to 1.0)
    """
    from numpy import array

    A = filter_a.graph.get_raw_tables()
    B = filter_b.graph.get_raw_tables()
    distance = 0
    for q, p in zip(A, B):
        a = array(q, copy=False)
        b = array(p, copy=False)
        for i in map(lambda x: randint(0, len(a)), range(n)):
            distance += sum(map(int,
                            [not bool((a[i] >> j) & 1) ^ bool((b[i] >> j) & 1)
                             for j in range(8)]))
    return distance / (8.0 * len(A) * n)


def convert_cmd(name, backend):
    "Convert an SBT to use a different back end."
    from .sbtmh import SigLeaf

    options = backend.split('(')
    backend = options.pop(0)
    backend = backend.lower().strip("'")

    if options:
      print(options)
      options = options[0].split(')')
      options = [options.pop(0)]
      #options = {}
    else:
      options = []

    if backend.lower() in ('ipfs', 'ipfsstorage'):
        backend = IPFSStorage
    elif backend.lower() in ('redis', 'redisstorage'):
        backend = RedisStorage
    elif backend.lower() in ('zip', 'zipstorage'):
        backend = ZipStorage
    elif backend.lower() in ('fs', 'fsstorage'):
        backend = FSStorage
        if options:
            options = [os.path.dirname(options[0]), os.path.basename(options[0])]
        else:
            # this is the default for SBT v2
            tag = '.sbt.' + os.path.basename(name)
            if tag.endswith('.sbt.json'):
                tag = tag[:-9]
            path = os.path.dirname(name)
            options = [path, tag]

    else:
        error('backend not recognized: {}'.format(backend))

    with backend(*options) as storage:
        sbt = SBT.load(name, leaf_loader=SigLeaf.load)
        sbt.save(name, storage=storage)
