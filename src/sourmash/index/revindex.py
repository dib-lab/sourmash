import weakref

from sourmash.index import Index, IndexSearchResult
from sourmash.minhash import MinHash
from sourmash.signature import SourmashSignature
from sourmash._lowlevel import ffi, lib
from sourmash.utils import RustObject, rustcall, decode_str, encode_str


class RevIndex(RustObject, Index):
    __dealloc_func__ = lib.revindex_free

    def __init__(
        self,
        *,
        signatures=None,
        signature_paths=None,
        template=None,
        threshold=0,
        queries=None,
        keep_sigs=False,
    ):
        self.template = template
        self.threshold = threshold
        self.queries = queries
        self.keep_sigs = keep_sigs
        self.signature_paths = signature_paths
        self._signatures = signatures

        if signature_paths is None or signatures is None:
            # delay initialization
            self._objptr = ffi.NULL
        else:
            self._init_inner()

    def _init_inner(self):
        if self._objptr != ffi.NULL:
            # Already initialized
            return

        if (
            self.signature_paths is None
            and not self._signatures
            and self._objptr == ffi.NULL
        ):
            raise ValueError("No signatures provided")
        elif (self.signature_paths or self._signatures) and self._objptr != ffi.NULL:
            raise NotImplementedError("Need to update RevIndex")

        attached_refs = weakref.WeakKeyDictionary()

        queries_ptr = ffi.NULL
        queries_size = 0
        if self.queries:
            # get list of rust objects
            collected = []
            for obj in queries:
                rv = obj._get_objptr()
                attached_refs[rv] = obj
                collected.append(rv)
            queries_ptr = ffi.new("SourmashSignature*[]", collected)
            queries_size = len(queries)

        template_ptr = ffi.NULL
        if self.template:
            if isinstance(self.template, MinHash):
                template_ptr = self.template._get_objptr()
            else:
                raise ValueError("Template must be a MinHash")

        search_sigs_ptr = ffi.NULL
        sigs_size = 0
        collected = []
        if self.signature_paths:
            for path in self.signature_paths:
                collected.append(encode_str(path))
            search_sigs_ptr = ffi.new("SourmashStr*[]", collected)
            sigs_size = len(signature_paths)

            self._objptr = rustcall(
                lib.revindex_new_with_paths,
                search_sigs_ptr,
                sigs_size,
                template_ptr,
                self.threshold,
                queries_ptr,
                queries_size,
                self.keep_sigs,
            )
        elif self._signatures:
            # force keep_sigs=True, and pass SourmashSignature directly to RevIndex.
            for sig in self._signatures:
                collected.append(sig._get_objptr())
            search_sigs_ptr = ffi.new("SourmashSignature*[]", collected)
            sigs_size = len(self._signatures)

            self._objptr = rustcall(
                lib.revindex_new_with_sigs,
                search_sigs_ptr,
                sigs_size,
                template_ptr,
                self.threshold,
                queries_ptr,
                queries_size,
            )

    def signatures(self):
        self._init_inner()

        size = ffi.new("uintptr_t *")
        sigs_ptr = self._methodcall(lib.revindex_signatures, size)
        size = size[0]

        sigs = []
        for i in range(size):
            sig = SourmashSignature._from_objptr(sigs_ptr[i])
            sigs.append(sig)

        for sig in sigs:
            yield sig

        #if self._signatures:
        #    yield from self._signatures
        #else:
        #    raise NotImplementedError("Call into Rust and retrieve sigs")

    def __len__(self):
        pass

    def insert(self, node):
        if self._signatures is None:
            self._signatures = []
        self._signatures.append(node)

    def save(self, path):
        pass

    @classmethod
    def load(cls, location):
        pass

    def select(self, ksize=None, moltype=None):
        if self.template:
            if ksize:
                self.template.ksize = ksize
            if moltype:
                self.template.moltype = moltype
        else:
            # TODO: deal with None/default values
            self.template = MinHash(ksize=ksize, moltype=moltype)

#    def search(self, query, *args, **kwargs):
#        """Return set of matches with similarity above 'threshold'.
#
#        Results will be sorted by similarity, highest to lowest.
#
#        Optional arguments:
#          * do_containment: default False. If True, use Jaccard containment.
#          * ignore_abundance: default False. If True, and query signature
#            and database support k-mer abundances, ignore those abundances.
#
#        Note, the "best only" hint is ignored by LCA_Database
#        """
#        if not query.minhash:
#            return []
#
#        # check arguments
#        if "threshold" not in kwargs:
#            raise TypeError("'search' requires 'threshold'")
#        threshold = kwargs["threshold"]
#        do_containment = kwargs.get("do_containment", False)
#        ignore_abundance = kwargs.get("ignore_abundance", False)
#
#        self._init_inner()
#
#        size = ffi.new("uintptr_t *")
#        results_ptr = self._methodcall(
#            lib.revindex_search,
#            query._get_objptr(),
#            threshold,
#            do_containment,
#            ignore_abundance,
#            size,
#        )
#
#        size = size[0]
#        if size == 0:
#            return []
#
#        results = []
#        for i in range(size):
#            match = SearchResult._from_objptr(results_ptr[i])
#            if match.score >= threshold:
#                results.append(IndexSearchResult(match.score, match.signature, match.filename))
#
#        return results
#
#    def gather(self, query, *args, **kwargs):
#        "Return the match with the best Jaccard containment in the database."
#        if not query.minhash:
#            return []
#
#        self._init_inner()
#
#        threshold_bp = kwargs.get("threshold_bp", 0.0)
#        threshold = threshold_bp / (len(query.minhash) * self.scaled)
#
#        results = []
#        size = ffi.new("uintptr_t *")
#        results_ptr = self._methodcall(
#            lib.revindex_gather, query._get_objptr(), threshold, True, True, size
#        )
#        size = size[0]
#        if size == 0:
#            return []
#
#        results = []
#        for i in range(size):
#            match = SearchResult._from_objptr(results_ptr[i])
#            if match.score >= threshold:
#                results.append(IndexSearchResult(match.score, match.signature, match.filename))
#
#        results.sort(reverse=True,
#                     key=lambda x: (x.score, x.signature.md5sum()))
#
#        return results[:1]

    @property
    def scaled(self):
        return self._methodcall(lib.revindex_scaled)


class SearchResult(RustObject):
    __dealloc_func__ = lib.searchresult_free

    @property
    def score(self):
        return self._methodcall(lib.searchresult_score)

    @property
    def signature(self):
        sig_ptr = self._methodcall(lib.searchresult_signature)
        return SourmashSignature._from_objptr(sig_ptr)

    @property
    def filename(self):
        result = decode_str(self._methodcall(lib.searchresult_filename))
        if result == "":
            return None
        return result
