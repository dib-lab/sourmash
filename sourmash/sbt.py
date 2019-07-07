#!/usr/bin/env python
"""
An implementation of sequence bloom trees, Solomon & Kingsford, 2015.

To try it out, do::

    factory = GraphFactory(ksize, tablesizes, n_tables)
    root = Node(factory)

    graph1 = factory()
    # ... add stuff to graph1 ...
    leaf1 = Leaf("a", graph1)
    root.add_node(leaf1)

For example, ::

    # filenames: list of fa/fq files
    # ksize: k-mer size
    # tablesizes: Bloom filter table sizes
    # n_tables: Number of tables

    factory = GraphFactory(ksize, tablesizes, n_tables)
    root = Node(factory)

    for filename in filenames:
        graph = factory()
        graph.consume_fasta(filename)
        leaf = Leaf(filename, graph)
        root.add_node(leaf)

then define a search function, ::

    def kmers(k, seq):
        for start in range(len(seq) - k + 1):
            yield seq[start:start + k]

    def search_transcript(node, seq, threshold):
        presence = [ node.data.get(kmer) for kmer in kmers(ksize, seq) ]
        if sum(presence) >= int(threshold * len(seq)):
            return 1
        return 0
"""

from __future__ import print_function, unicode_literals, division

from collections import namedtuple
try:
    from collections.abc import Mapping
except ImportError:  # Python 2...
    from collections import Mapping

from copy import copy
import json
import math
import os
from random import randint, random
import sys
from tempfile import NamedTemporaryFile
import zipfile

import khmer

try:
    load_nodegraph = khmer.load_nodegraph
except AttributeError:
    load_nodegraph = khmer.Nodegraph.load

from .sbt_storage import FSStorage, TarStorage, IPFSStorage, RedisStorage, ZipStorage
from .logging import error, notify, debug


STORAGES = {
    'TarStorage': TarStorage,
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
        return khmer.Nodegraph(self.ksize, self.starting_size, self.n_tables)

    def init_args(self):
        return (self.ksize, self.starting_size, self.n_tables)


class SBT(object):
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

    Notes
    -----
    We use two dicts to store the tree structure: One for the internal nodes,
    and another for the leaves (datasets).
    """

    def __init__(self, factory, d=2, storage=None):
        self.factory = factory
        self._nodes = {}
        self._missing_nodes = set()
        self._leaves = {}
        self.d = d
        self.next_node = 0
        self.storage = storage
        self.is_ready = False

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

    def add_node(self, leaf, update_internal=True):
        pos = self.new_node_pos(leaf)

        if pos == 0:  # empty tree; initialize w/node.
            n = Node(self.factory, name="internal." + str(pos))
            self._nodes[0] = n
            pos = self.new_node_pos(leaf)

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
            self._leaves[c2.pos] = leaf
            del self._leaves[p.pos]

            if update_internal:
                for child in (p.node, leaf):
                    child.update(n)
            else:
               self.is_ready = False
        elif isinstance(p.node, Node):
            self._leaves[pos] = leaf
            if update_internal:
                leaf.update(p.node)
            else:
               self.is_ready = False
        elif p.node is None:
            n = Node(self.factory, name="internal." + str(p.pos))
            self._nodes[p.pos] = n
            c1 = self.children(p.pos)[0]
            self._leaves[c1.pos] = leaf
            if update_internal:
                leaf.update(n)
            else:
               self.is_ready = False

        if update_internal:
            # update all parents!
            p = self.parent(p.pos)
            while p:
                self._rebuild_node(p.pos)
                leaf.update(self._nodes[p.pos])
                p = self.parent(p.pos)
        else:
           self.is_ready = False

    def find(self, search_fn, *args, **kwargs):
        "Search the tree using `search_fn`."

        # initialize search queue with top node of tree
        if not self.is_ready:
            self._fill_internal()

        matches = []
        visited, queue = set(), [0]

        # while the queue is not empty, load each node and apply search
        # function.
        while queue:
            node_p = queue.pop(0)

            # repair while searching.
            node_g = self._leaves.get(node_p, None)
            if node_g is None:
                node_g = self._nodes.get(node_p, None)

            if node_g is None:
                if node_p in self._missing_nodes:
                    self._rebuild_node(node_p)
                    node_g = self._nodes[node_p]
                else:
                    continue

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
        return matches

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
        version = 5

        if path.endswith('.sbt.json'):
            path = path[:-9]
        fn = os.path.abspath(path + '.sbt.json')

        if storage is None:
            # default storage
            location = os.path.dirname(fn)
            subdir = '.sbt.{}'.format(os.path.basename(path))

            storage = FSStorage(location, subdir)
            fn = os.path.join(location, fn)

        backend = [k for (k, v) in STORAGES.items() if v == type(storage)][0]

        info = {}
        info['d'] = self.d
        info['version'] = version
        info['storage'] = {
            'backend': backend,
            'args': storage.init_args()
        }
        info['factory'] = {
            'class': GraphFactory.__name__,
            'args': self.factory.init_args()
        }

        if not self.is_ready and structure_only is False:
            self._fill_internal_and_save(storage, sparseness)

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

                data['filename'] = node.save(data['filename'])

            if isinstance(node, Node):
                nodes[i] = data
            else:
                leaves[i] = data

            notify("{} of {} nodes saved".format(n+1, total_nodes), end='\r')

        notify("\nFinished saving nodes, now saving SBT json file.")
        info['nodes'] = nodes
        info['leaves'] = leaves
        with open(fn, 'w') as fp:
            json.dump(info, fp)

        return fn

    @classmethod
    def load(cls, location, leaf_loader=None, storage=None, print_version_warning=True):
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
        if zipfile.is_zipfile(location):
            tempfile = NamedTemporaryFile()
            with zipfile.ZipFile(location, 'r') as zf:
                tempfile.write(zf.read('tree.sbt.json'))
                tempfile.flush()

            dirname = os.path.dirname(tempfile.name)
            sbt_name = os.path.basename(tempfile.name)
            storage = ZipStorage(location)
        else:
            dirname = os.path.dirname(os.path.abspath(location))
            sbt_name = os.path.basename(location)
            if sbt_name.endswith('.sbt.json'):
                sbt_name = sbt_name[:-9]

        loaders = {
            1: cls._load_v1,
            2: cls._load_v2,
            3: cls._load_v3,
            4: cls._load_v4,
            5: cls._load_v5,
        }

        # @CTB hack: check to make sure khmer Nodegraph supports the
        # correct methods.
        x = khmer.Nodegraph(1, 1, 1)
        try:
            x.count(10)
        except TypeError:
            raise Exception("khmer version is too old; need >= 2.1,<3")

        if leaf_loader is None:
            leaf_loader = Leaf.load

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

        if version < 3 and storage is None:
            storage = FSStorage(dirname, '.sbt.{}'.format(sbt_name))

        return loaders[version](jnodes, leaf_loader, dirname, storage,
                                print_version_warning)

    @staticmethod
    def _load_v1(jnodes, leaf_loader, dirname, storage, print_version_warning=True):

        if jnodes[0] is None:
            raise ValueError("Empty tree!")

        sbt_nodes = {}
        sbt_leaves = {}

        max_node = 0

        sample_bf = os.path.join(dirname, jnodes[0]['filename'])
        ksize, tablesize, ntables = khmer.extract_nodegraph_info(sample_bf)[:3]
        factory = GraphFactory(ksize, tablesize, ntables)

        for k, jnode in enumerate(jnodes):
            if jnode is None:
                continue

            jnode['filename'] = os.path.join(dirname, jnode['filename'])

            if 'internal' in jnode['name']:
                jnode['factory'] = factory
                sbt_node = Node.load(jnode, storage)
                sbt_nodes[k] = sbt_node
            else:
                sbt_node = leaf_loader(jnode, storage)
                sbt_leaves[k] = sbt_node

            max_node = max(max_node, k)

        tree = SBT(factory)
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
    def _load_v2(cls, info, leaf_loader, dirname, storage, print_version_warning=True):
        nodes = {int(k): v for (k, v) in info['nodes'].items()}

        if nodes[0] is None:
            raise ValueError("Empty tree!")

        sbt_nodes = {}
        sbt_leaves = {}

        max_node = 0

        sample_bf = os.path.join(dirname, nodes[0]['filename'])
        k, size, ntables = khmer.extract_nodegraph_info(sample_bf)[:3]
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

            max_node = max(max_node, k)

        tree = cls(factory, d=info['d'])
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
    def _load_v3(cls, info, leaf_loader, dirname, storage, print_version_warning=True):
        nodes = {int(k): v for (k, v) in info['nodes'].items()}

        if not nodes:
            raise ValueError("Empty tree!")

        sbt_nodes = {}
        sbt_leaves = {}

        klass = STORAGES[info['storage']['backend']]
        if info['storage']['backend'] == "FSStorage":
            storage = FSStorage(dirname, info['storage']['args']['path'])
        elif storage is None:
            storage = klass(**info['storage']['args'])

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

        tree = cls(factory, d=info['d'], storage=storage)
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
    def _load_v4(cls, info, leaf_loader, dirname, storage, print_version_warning=True):
        nodes = {int(k): v for (k, v) in info['nodes'].items()}

        if not nodes:
            raise ValueError("Empty tree!")

        sbt_nodes = {}
        sbt_leaves = {}

        klass = STORAGES[info['storage']['backend']]
        if info['storage']['backend'] == "FSStorage":
            storage = FSStorage(dirname, info['storage']['args']['path'])
        elif storage is None:
            storage = klass(**info['storage']['args'])

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

        tree = cls(factory, d=info['d'], storage=storage)
        tree._nodes = sbt_nodes
        tree._leaves = sbt_leaves
        tree._missing_nodes = {i for i in range(max_node)
                               if i not in sbt_nodes and i not in sbt_leaves}

        if print_version_warning:
            error("WARNING: this is an old index version, please run `sourmash migrate` to update it.")
            error("WARNING: proceeding with execution, but it will take longer to finish!")

        return tree

    @classmethod
    def _load_v5(cls, info, leaf_loader, dirname, storage, print_version_warning=True):
        nodes = {}
        if 'nodes' in info:
            nodes = {int(k): v for (k, v) in info['nodes'].items()}

        leaves = {int(k): v for (k, v) in info['leaves'].items()}

        if not leaves:
            raise ValueError("Empty tree!")

        sbt_nodes = {}
        sbt_leaves = {}

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

        tree = cls(factory, d=info['d'], storage=storage)
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

    def _fill_internal_and_save(self, storage, sparseness=0.0):

        def fill_nodegraphs_and_save(node, *args, **kwargs):
            children = kwargs['children']
            for child in children:
                if child.node is not None:
                    child.node.update(node)

                    if isinstance(node, Node) and random() - sparseness > 0:
                        child.node.storage = storage
                        child.node.save(os.path.basename(node.name))

                        child.node.unload()
            return True

        self._fill_up(fill_nodegraphs_and_save)
        self.is_ready = True

    def _fill_internal(self):

        def fill_nodegraphs(node, *args, **kwargs):
            children = kwargs['children']
            for child in children:
                if child.node is not None:
                    child.node.update(node)
            return True

        self._fill_up(fill_nodegraphs)
        self.is_ready = True

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
                fpr=khmer.calc_expected_collisions(self.data, True, 1.1))

    def save(self, path):
        # We need to do this tempfile dance because khmer only load
        # data from files.
        with NamedTemporaryFile(suffix=".gz") as f:
            self.data.save(f.name)
            f.file.flush()
            f.file.seek(0)
            return self.storage.save(path, f.read())

    @property
    def data(self):
        if self._data is None:
            if self._path is None:
                self._data = self._factory()
            else:
                data = self.storage.load(self._path)
                # We need to do this tempfile dance because khmer only load
                # data from files.
                with NamedTemporaryFile(suffix=".gz") as f:
                    f.write(data)
                    f.file.flush()
                    self._data = load_nodegraph(f.name)
        return self._data

    @data.setter
    def data(self, new_data):
        self._data = new_data

    @staticmethod
    def load(info, storage=None):
        new_node = Node(info['factory'],
                        name=info['name'],
                        path=info['filename'],
                        storage=storage)
        new_node.metadata = info.get('metadata', {})
        return new_node

    def unload(self):
        pass

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
                fpr=khmer.calc_expected_collisions(self.data, True, 1.1))

    @property
    def data(self):
        if self._data is None:
            data = self.storage.load(self._path)
            # We need to do this tempfile dance because khmer only load
            # data from files.
            with NamedTemporaryFile(suffix=".gz") as f:
                f.write(data)
                f.file.flush()
                self._data = load_nodegraph(f.name)
        return self._data

    @data.setter
    def data(self, new_data):
        self._data = new_data

    def save(self, path):
        # We need to do this tempfile dance because khmer only load
        # data from files.
        with NamedTemporaryFile(suffix=".gz") as f:
            self.data.save(f.name)
            f.file.flush()
            f.file.seek(0)
            return self.storage.save(path, f.read())

    def update(self, parent):
        parent.data.update(self.data)

    @classmethod
    def load(cls, info, storage=None):
        return cls(info['metadata'],
                   name=info['name'],
                   path=info['filename'],
                   storage=storage)

    def unload(self):
        pass


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


def parse_backend_args(name, backend):
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
    elif backend.lower() in ('tar', 'tarstorage'):
        backend = TarStorage
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

    return backend, options


def convert_cmd(name, backend_args):
    from .sbtmh import SigLeaf

    backend, options = parse_backend_args(name, backend_args)

    with backend(*options) as storage:
        sbt = SBT.load(name, leaf_loader=SigLeaf.load)
        sbt.save(name, storage=storage)
