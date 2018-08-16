from __future__ import print_function
from __future__ import division

import collections
from io import BytesIO, TextIOWrapper
import sys

from .sbt import Leaf, SBT, GraphFactory
from . import signature


def load_sbt_index(filename, print_version_warning=True):
    "Load and return an SBT index."
    return SBT.load(filename, leaf_loader=SigLeaf.load,
                    print_version_warning=print_version_warning)


def create_sbt_index(bloom_filter_size=1e5, n_children=2, ksize=1):
    "Create an empty SBT index."
    factory = GraphFactory(ksize, bloom_filter_size, 4)
    tree = SBT(factory, d=n_children)
    return tree


def search_sbt_index(tree, query, threshold):
    """\
    Search an SBT index `tree` with signature `query` for matches above
    `threshold`.

    Usage:

        for match_sig, similarity in search_sbt_index(tree, query, threshold):
           ...
    """
    for leaf in tree.find(search_minhashes, query, threshold):
        similarity = query.similarity(leaf.data)
        yield leaf.data, similarity


class SigLeaf(Leaf):
    def __str__(self):
        return '**Leaf:{name} -> {metadata}'.format(
                name=self.name, metadata=self.metadata)

    def save(self, path):
        # this is here only for triggering the property load
        # before we reopen the file (and overwrite the previous
        # content...)
        self.data

        buf = BytesIO()
        with TextIOWrapper(buf) as out:
            sigs = self.data
            if not isinstance(self.data, collections.Sequence):
                sigs = [self.data]

            signature.save_signatures(sigs, out)
            out.flush()
            return self.storage.save(path, buf.getvalue())

    def update(self, parent):
        if isinstance(self.data, collections.Sequence):
            # TODO: select proper signature
            minhash = self.data[0].minhash
        else:
            minhash = self.data.minhash

        for v in minhash.get_mins():
            parent.data.count(v)
        min_n_below = parent.metadata.get('min_n_below', sys.maxsize)
        min_n_below = min(len(minhash), min_n_below)

        if min_n_below == 0:
            min_n_below = 1

        parent.metadata['min_n_below'] = min_n_below

    @property
    def data(self):
        if self._data is None:
            buf = BytesIO(self.storage.load(self._path))
            with TextIOWrapper(buf) as data:
                self._data = list(signature.load_signatures(data))
        return self._data

    @data.setter
    def data(self, new_data):
        self._data = new_data


### Search functionality.

def _max_jaccard_underneath_internal_node(node, hashes):
    """\
    calculate the maximum possibility similarity score below
    this node, based on the number of matches in 'hashes' at this node,
    divided by the smallest minhash size below this node.

    This should yield be an upper bound on the Jaccard similarity
    for any signature below this point.
    """
    if len(hashes) == 0:
        return 0.0

    # count the maximum number of hash matches beneath this node
    get = node.data.get
    matches = sum(1 for value in hashes if get(value))

    # get the size of the smallest collection of hashes below this point
    min_n_below = node.metadata.get('min_n_below', -1)

    if min_n_below == -1:
        raise Exception('cannot do similarity search on this SBT; need to rebuild.')

    # max of numerator divided by min of denominator => max Jaccard
    max_score = float(matches) / min_n_below

    return max_score


def select_signature(node, sig):
    if isinstance(node.data, collections.Sequence):
        # TODO: select proper signature
        matches = [q for q in node.data if q.minhash.ksize == sig.minhash.ksize]
        if matches:
            return matches[0]
        return None
    return node.data


def search_minhashes(node, sig, threshold, results=None, downsample=True):
    """\
    Default tree search function, searching for best Jaccard similarity.
    """
    mins = sig.minhash.get_mins()
    score = 0

    if isinstance(node, SigLeaf):
        node_mh = select_signature(node, sig).minhash
        try:
            score = node_mh.similarity(sig.minhash)
        except Exception as e:
            if 'mismatch in max_hash' in str(e) and downsample:
                xx = sig.minhash.downsample_max_hash(node_mh)
                yy = node_mh.downsample_max_hash(sig.minhash)

                score = yy.similarity(xx)
            else:
                raise

    else:  # Node minhash comparison
        score = _max_jaccard_underneath_internal_node(node, mins)

    if results is not None:
        results[node.name] = score

    if score >= threshold:
        return 1

    return 0


class SearchMinHashesFindBest(object):
    def __init__(self, downsample=True):
        self.best_match = 0.
        self.downsample = downsample

    def search(self, node, sig, threshold, results=None):
        mins = sig.minhash.get_mins()
        score = 0

        if isinstance(node, SigLeaf):
            node_mh = select_signature(node, sig).minhash
            try:
                score = node_mh.similarity(sig.minhash)
            except Exception as e:
                if 'mismatch in max_hash' in str(e) and self.downsample:
                    xx = sig.minhash.downsample_max_hash(node_mh)
                    yy = node_mh.downsample_max_hash(sig.minhash)

                    score = yy.similarity(xx)
                else:
                    raise
        else:  # internal object, not leaf.
            score = _max_jaccard_underneath_internal_node(node, mins)

        if results is not None:
            results[node.name] = score

        if score >= threshold:
            # have we done better than this elsewhere? if yes, truncate.
            if score > self.best_match:
                # update best if it's a leaf node...
                if isinstance(node, SigLeaf):
                    self.best_match = score
                return 1

        return 0


def search_minhashes_containment(node, sig, threshold,
                                 results=None, downsample=True):
    mins = sig.minhash.get_mins()

    if isinstance(node, SigLeaf):
        node_mh = select_signature(node, sig).minhash
        try:
            matches = node_mh.count_common(sig.minhash)
        except Exception as e:
            if 'mismatch in max_hash' in str(e) and downsample:
                xx = sig.minhash.downsample_max_hash(node_mh)
                yy = node_mh.downsample_max_hash(sig.minhash)

                matches = yy.count_common(xx)
            else:
                raise

    else:  # Node or Leaf, Nodegraph by minhash comparison
        get = node.data.get
        matches = sum(1 for value in mins if get(value))

    if results is not None:
        results[node.name] = float(matches) / len(mins)

    if len(mins) and float(matches) / len(mins) >= threshold:
        return 1
    return 0


class SearchMinHashesFindBestIgnoreMaxHash(object):
    def __init__(self):
        self.best_match = 0.

    def search(self, node, sig, threshold, results=None):
        mins = sig.minhash.get_mins()

        if isinstance(node, SigLeaf):
            node_mh = select_signature(node, sig).minhash
            max_scaled = max(node_mh.scaled, sig.minhash.scaled)

            mh1 = node_mh.downsample_scaled(max_scaled)
            mh2 = sig.minhash.downsample_scaled(max_scaled)
            matches = mh1.count_common(mh2)
        else:  # Node or Leaf, Nodegraph by minhash comparison
            get = node.data.get
            matches = sum(1 for value in mins if get(value))

        score = 0
        if not len(mins):
            return 0

        score = float(matches) / len(mins)

        if results is not None:
            results[node.name] = score

        if score >= threshold:
            # have we done better than this? if yes, truncate.
            if float(matches) / len(mins) > self.best_match:
                # update best if it's a leaf node...
                if isinstance(node, SigLeaf):
                    self.best_match = float(matches) / len(mins)
                return 1

        return 0
