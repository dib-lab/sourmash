from io import BytesIO
import sys

from .sbt import Leaf, SBT, GraphFactory
from . import signature


def load_sbt_index(filename, *, print_version_warning=True, cache_size=None):
    "Load and return an SBT index."
    return SBT.load(filename, leaf_loader=SigLeaf.load,
                    print_version_warning=print_version_warning,
                    cache_size=cache_size)


def create_sbt_index(bloom_filter_size=1e5, n_children=2):
    "Create an empty SBT index."
    factory = GraphFactory(1, bloom_filter_size, 4)
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
    for leaf in tree.find(search_minhashes, query, threshold, unload_data=True):
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

        buf = signature.save_signatures([self.data], compression=1)
        return self.storage.save(path, buf)

    def update(self, parent):
        mh = self.data.minhash
        parent.data.update(mh)

    @property
    def data(self):
        if self._data is None:
            buf = BytesIO(self.storage.load(self._path))
            self._data = signature.load_one_signature(buf)
        return self._data

    @data.setter
    def data(self, new_data):
        self._data = new_data


### Search functionality.

def _max_jaccard_underneath_internal_node(node, query):
    """\
    calculate the maximum possibility similarity score below
    this node, based on the number of matches in 'hashes' at this node,
    divided by the size of the query.

    This should yield be an upper bound on the Jaccard similarity
    for any signature below this point.
    """
    mh = query.minhash

    if len(mh) == 0:
        return 0.0

    if mh.track_abundance:
        # In this case we need to use the upper bound for angular similarity
        max_score = node.data.angular_similarity_upper_bound(mh)
    else:
        # In this case we are working with similarity/containment:
        # J(A, B) = |A intersection B| / |A union B|
        # If we use only |A| as denominator, it is the containment
        # Because |A| <= |A union B|, it is also an upper bound on the max jaccard

        # count the maximum number of hash matches beneath this node
        matches = node.data.matches(mh)

        max_score = float(matches) / len(mh)

    return max_score


def search_minhashes(node, sig, threshold, results=None):
    """\
    Default tree search function, searching for best Jaccard similarity.
    """
    assert results is None

    score = 0

    if isinstance(node, SigLeaf):
        score = node.data.minhash.similarity(sig.minhash)
    else:  # Node minhash comparison
        score = _max_jaccard_underneath_internal_node(node, sig)

    if score >= threshold:
        return 1

    return 0


class SearchMinHashesFindBest(object):
    def __init__(self):
        self.best_match = 0.

    def search(self, node, sig, threshold, results=None):
        assert results is None
        score = 0

        if isinstance(node, SigLeaf):
            score = node.data.minhash.similarity(sig.minhash)
        else:  # internal object, not leaf.
            score = _max_jaccard_underneath_internal_node(node, sig)

        if score >= threshold:
            # have we done better than this elsewhere? if yes, truncate.
            if score > self.best_match:
                # update best if it's a leaf node...
                if isinstance(node, SigLeaf):
                    self.best_match = score
                return 1

        return 0


def search_minhashes_containment(node, sig, threshold, results=None, downsample=True):
    assert results is None
    mh = sig.minhash

    if isinstance(node, SigLeaf):
        matches = node.data.minhash.count_common(mh, downsample)
    else:  # Node or Leaf, Nodegraph by minhash comparison
        matches = node.data.matches(mh)

    if len(mh) and float(matches) / len(mh) >= threshold:
        return 1
    return 0


def search_minhashes_max_containment(node, sig, threshold, results=None,
                                     downsample=True):
    assert results is None

    mh = sig.minhash

    if isinstance(node, SigLeaf):
        node_mh = node.data.minhash

        matches = node_mh.count_common(mh, downsample)
        node_size = len(node_mh)
    else:  # Node or Leaf, Nodegraph by minhash comparison
        matches = node.data.matches(mh)

        # get the size of the smallest collection of hashes below this point
        node_size = node.metadata.get('min_n_below', -1)

        if node_size == -1:
            raise Exception('cannot do max_containment search on this SBT; need to rebuild.')

    denom = min((len(mh), node_size))

    if len(mh) and matches / denom >= threshold:
        return 1

    return 0


class GatherMinHashes(object):
    def __init__(self):
        self.best_match = 0

    def search(self, node, query, threshold, results=None):
        assert results is None

        mh = query.minhash
        if not len(mh):
            return 0

        if isinstance(node, SigLeaf):
            matches = mh.count_common(node.data.minhash, True)
        else:  # Nodegraph by minhash comparison
            matches = node.data.matches(mh)

        if not matches:
            return 0

        score = float(matches) / len(mh)

        if score < threshold:
            return 0

        # have we done better than this? if no, truncate searches below.
        if score >= self.best_match:
            # update best if it's a leaf node...
            if isinstance(node, SigLeaf):
                self.best_match = score
            return 1

        return 0
