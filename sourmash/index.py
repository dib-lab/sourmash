"An Abstract Base Class for collections of signatures."

from abc import ABCMeta, abstractmethod
from collections import namedtuple

# compatible with Python 2 *and* 3:
ABC = ABCMeta("ABC", (object,), {"__slots__": ()})


class Index(ABC):
    @abstractmethod
    def signatures(self):
        "Return an iterator over all signatures in the Index object."

    @abstractmethod
    def insert(self, signature):
        """ """

    @abstractmethod
    def save(self, path, storage=None, sparseness=0.0, structure_only=False):
        """ """

    @classmethod
    @abstractmethod
    def load(cls, location, leaf_loader=None, storage=None, print_version_warning=True):
        """ """

    def find(self, search_fn, *args, **kwargs):
        """Use search_fn to find matching signatures in the index.

        search_fn(other_sig, *args) should return a boolean that indicates
        whether other_sig is a match.

        Returns a list.
        """

        matches = []

        for node in self.signatures():
            if search_fn(node, *args):
                matches.append(node)
        return matches

    def search(self, query, *args, **kwargs):
        """Return set of matches with similarity above 'threshold'.

        Results will be sorted by similarity, highest to lowest.

        Optional arguments accepted by all Index subclasses:
          * do_containment: default False. If True, use Jaccard containment.
          * best_only: default False. If True, allow optimizations that
            may. May discard matches better than threshold, but first match
            is guaranteed to be best.
          * ignore_abundance: default False. If True, and query signature
            and database support k-mer abundances, ignore those abundances.

        Note, the "best only" hint is ignored by LinearIndex.
        """

        # check arguments
        if 'threshold' not in kwargs:
            raise TypeError("'search' requires 'threshold'")
        threshold = kwargs['threshold']

        do_containment = kwargs.get('do_containment', False)
        ignore_abundance = kwargs.get('ignore_abundance', False)

        # configure search - containment? ignore abundance?
        if do_containment:
            query_match = lambda x: query.contained_by(x, downsample=True)
        else:
            query_match = lambda x: query.similarity(
                x, downsample=True, ignore_abundance=ignore_abundance)

        # do the actual search:
        matches = []

        for ss in self.signatures():
            similarity = query_match(ss)
            if similarity >= threshold:
                matches.append((similarity, ss, self.filename))

        # sort!
        matches.sort(key=lambda x: -x[0])
        return matches

    def gather(self, query, *args, **kwargs):
        "Return the match with the best Jaccard containment in the Index."
        results = []
        for ss in self.signatures():
            cont = query.minhash.contained_by(ss.minhash, True)
            if cont:
                results.append((cont, ss, self.filename))

        results.sort(reverse=True, key=lambda x: (x[0], x[1].name()))

        return results

    @abstractmethod
    def select(self, ksize=None, moltype=None):
        ""

class LinearIndex(Index):
    def __init__(self, _signatures=None, filename=None):
        self._signatures = []
        if _signatures:
            self._signatures = list(_signatures)
        self.filename = filename

    def signatures(self):
        return iter(self._signatures)

    def __len__(self):
        return len(self._signatures)

    def insert(self, node):
        self._signatures.append(node)

    def save(self, path):
        from .signature import save_signatures
        with open(path, 'wt') as fp:
            save_signatures(self.signatures(), fp)

    @classmethod
    def load(cls, location):
        from .signature import load_signatures
        si = load_signatures(location)

        lidx = LinearIndex(si, filename=location)
        return lidx

    def select(self, ksize=None, moltype=None):
        def select_sigs(siglist, ksize, moltype):
            for ss in siglist:
#                print(self.filename, ss, ss.minhash.ksize, ss.minhash.moltype,
#                      ksize, moltype)
                if (ksize is None or ss.minhash.ksize == ksize) and \
                   (moltype is None or ss.minhash.moltype == moltype):
                   yield ss

        siglist=select_sigs(self._signatures, ksize, moltype)
        return LinearIndex(siglist, self.filename)
