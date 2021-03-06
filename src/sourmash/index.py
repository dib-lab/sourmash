"An Abstract Base Class for collections of signatures."

from abc import abstractmethod, ABC
from collections import namedtuple, Counter


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
        if not query.minhash:             # empty query? quit.
            return []

        scaled = query.minhash.scaled
        if not scaled:
            raise ValueError('gather requires scaled signatures')

        threshold_bp = kwargs.get('threshold_bp', 0.0)
        threshold = 0.0

        # are we setting a threshold?
        if threshold_bp:
            # if we have a threshold_bp of N, then that amounts to N/scaled
            # hashes:
            n_threshold_hashes = float(threshold_bp) / scaled

            # that then requires the following containment:
            threshold = n_threshold_hashes / len(query.minhash)

            # is it too high to ever match? if so, exit.
            if threshold > 1.0:
                return []

        # actually do search!
        results = []
        for ss in self.signatures():
            cont = query.minhash.contained_by(ss.minhash, True)
            if cont and cont >= threshold:
                results.append((cont, ss, self.filename))

        results.sort(reverse=True, key=lambda x: (x[0], x[1].md5sum()))

        return results

    def counter_gather(self, query, *args, **kwargs):
        "Perform compositional analysis of the query using the gather algorithm"
        if not query.minhash:             # empty query? quit.
            return []

        scaled = query.minhash.scaled
        if not scaled:
            raise ValueError('gather requires scaled signatures')

        threshold_bp = kwargs.get('threshold_bp', 0.0)
        threshold = 0.0
        n_threshold_hashes = 0

        # are we setting a threshold?
        if threshold_bp:
            # if we have a threshold_bp of N, then that amounts to N/scaled
            # hashes:
            n_threshold_hashes = float(threshold_bp) / scaled

            # that then requires the following containment:
            threshold = n_threshold_hashes / len(query.minhash)

            # is it too high to ever match? if so, exit.
            if threshold > 1.0:
                return []

        # Pre-loading signatures so we can index datasets
        signatures = list(self.signatures())

        # Process all datasets and create a Counter containing the size
        # of hashes in common between query and each signature
        counter = Counter()
        for (i, ss) in enumerate(signatures):
            counter[i] = query.minhash.count_common(ss.minhash, True)

        # Decompose query into matching signatures using a greedy approach (gather)
        results = []
        match_size = n_threshold_hashes
        while counter and match_size >= n_threshold_hashes:
            most_common = counter.most_common()
            dataset_id, size = most_common[0]
            if size >= n_threshold_hashes:
                match_size = size
            else:
                break

            match = signatures[dataset_id]
            del counter[dataset_id]
            cont = query.minhash.contained_by(match.minhash, True)
            if cont and cont >= threshold:
                results.append((cont, match, getattr(self, "filename", None)))

            # Prepare counter for finding the next match by decrementing
            # all hashes found in the current match in other datasets
            for (dataset_id, _) in most_common:
                counter[dataset_id] -= signatures[dataset_id].minhash.count_common(match.minhash, True)
                if counter[dataset_id] == 0:
                    del counter[dataset_id]

        results.sort(reverse=True, key=lambda x: (x[0], x[1].md5sum()))

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
        si = load_signatures(location, do_raise=True)

        lidx = LinearIndex(si, filename=location)
        return lidx

    def select(self, ksize=None, moltype=None):
        def select_sigs(siglist, ksize, moltype):
            for ss in siglist:
                if (ksize is None or ss.minhash.ksize == ksize) and \
                   (moltype is None or ss.minhash.moltype == moltype):
                   yield ss

        siglist=select_sigs(self._signatures, ksize, moltype)
        return LinearIndex(siglist, self.filename)
