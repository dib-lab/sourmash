"An Abstract Base Class for collections of signatures."

import os
import sourmash
from abc import abstractmethod, ABC
from collections import namedtuple, Counter
import zipfile
import copy

from .search import make_jaccard_search_query, make_gather_query

# generic return tuple for Index.search and Index.gather
IndexSearchResult = namedtuple('Result', 'score, signature, location')

class Index(ABC):
    is_database = False

    @property
    def location(self):
        "Return a resolvable location for this index, if possible."
        return None

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
    def load(cls, location, leaf_loader=None, storage=None,
             print_version_warning=True):
        """ """

    def find(self, search_fn, query, **kwargs):
        """Use search_fn to find matching signatures in the index.

        search_fn follows the protocol in JaccardSearch objects.

        Returns a list.
        """
        # first: is this query compatible with this search?
        search_fn.check_is_compatible(query)

        # ok! continue!

        # this set of signatures may be heterogenous in scaled/num values;
        # define some processing functions to downsample appropriately.
        query_mh = query.minhash
        assert not query_mh.track_abundance
        if query_mh.scaled:
            # make query and subject compatible w/scaled.
            query_scaled = query_mh.scaled

            def prepare_subject(subj_mh):
                assert subj_mh.scaled
                if subj_mh.track_abundance:
                    subj_mh = subj_mh.flatten()

                # downsample subject to highest scaled
                subj_scaled = subj_mh.scaled
                if subj_scaled < query_scaled:
                    return subj_mh.downsample(scaled=query_scaled)
                else:
                    return subj_mh

            def prepare_query(query_mh, subj_mh):
                assert subj_mh.scaled

                # downsample query to highest scaled
                subj_scaled = subj_mh.scaled
                if subj_scaled > query_scaled:
                    return query_mh.downsample(scaled=subj_scaled)
                else:
                    return query_mh

        else:                   # num
            query_num = query_mh.num

            def prepare_subject(subj_mh):
                assert subj_mh.num
                if subj_mh.track_abundance:
                    subj_mh = subj_mh.flatten()

                # downsample subject to smallest num
                subj_num = subj_mh.num
                if subj_num > query_num:
                    return subj_mh.downsample(num=query_num)
                else:
                    return subj_mh

            def prepare_query(query_mh, subj_mh):
                assert subj_mh.num
                # downsample query to smallest num
                subj_num = subj_mh.num
                if subj_num < query_num:
                    return query_mh.downsample(num=subj_num)
                else:
                    return query_mh

        # now, do the search!
        for subj in self.signatures():
            subj_mh = prepare_subject(subj.minhash)
            # note: we run prepare_query here on the original query minhash.
            query_mh = prepare_query(query.minhash, subj_mh)

            assert not query_mh.track_abundance
            assert not subj_mh.track_abundance

            shared_size, total_size = query_mh.intersection_and_union_size(subj_mh)

            query_size = len(query_mh)
            subj_size = len(subj_mh)

            score = search_fn.score_fn(query_size,
                                       shared_size,
                                       subj_size,
                                       total_size)

            if search_fn.passes(score):
                # note: here we yield the original signature, not the
                # downsampled minhash.
                if search_fn.collect(score, subj):
                    yield subj, score

    def search_abund(self, query, *, threshold=None, **kwargs):
        """Return set of matches with angular similarity above 'threshold'.

        Results will be sorted by similarity, highest to lowest.
        """
        if not query.minhash.track_abundance:
            raise TypeError("'search_abund' requires query signature with abundance information")

        # check arguments
        if threshold is None:
            raise TypeError("'search_abund' requires 'threshold'")
        threshold = float(threshold)

        # do the actual search:
        matches = []
        for subj in self.signatures():
            if not subj.minhash.track_abundance:
                raise TypeError("'search_abund' requires subject signatures with abundance information")
            score = query.similarity(subj)
            if score >= threshold:
                matches.append(IndexSearchResult(score, subj, self.location))

        # sort!
        matches.sort(key=lambda x: -x.score)
        return matches

    def search(self, query, *, threshold=None,
               do_containment=False, do_max_containment=False,
               best_only=False, **kwargs):
        """Return set of matches with similarity above 'threshold'.

        Results will be sorted by similarity, highest to lowest.

        Optional arguments accepted by all Index subclasses:
          * do_containment: default False. If True, use Jaccard containment.
          * best_only: default False. If True, allow optimizations that
            may. May discard matches better than threshold, but first match
            is guaranteed to be best.
        """
        # check arguments
        if threshold is None:
            raise TypeError("'search' requires 'threshold'")
        threshold = float(threshold)

        search_obj = make_jaccard_search_query(do_containment=do_containment,
                                               do_max_containment=do_max_containment,
                                               best_only=best_only,
                                               threshold=threshold)

        # do the actual search:
        matches = []

        for subj, score in self.find(search_obj, query, **kwargs):
            matches.append(IndexSearchResult(score, subj, self.location))

        # sort!
        matches.sort(key=lambda x: -x.score)
        return matches

    def prefetch(self, query, threshold_bp, **kwargs):
        "Return all matches with minimum overlap."
        query_mh = query.minhash
        scaled = query_mh.scaled

        if not self:            # empty database? quit.
            raise ValueError("no signatures to search")

        search_fn = make_gather_query(query.minhash, threshold_bp,
                                      best_only=False)

        for subj, score in self.find(search_fn, query, **kwargs):
            yield IndexSearchResult(score, subj, self.location)

    def gather(self, query, threshold_bp=None, **kwargs):
        "Return the match with the best Jaccard containment in the Index."

        results = []
        for result in self.prefetch(query, threshold_bp, **kwargs):
            results.append(result)

        # sort results by best score.
        results.sort(reverse=True,
                     key=lambda x: (x.score, x.signature.md5sum()))

        return results[:1]

    @abstractmethod
    def select(self, ksize=None, moltype=None, scaled=None, num=None,
               abund=None, containment=None):
        """Return Index containing only signatures that match requirements.

        Current arguments can be any or all of:
        * ksize
        * moltype
        * scaled
        * num
        * containment

        'select' will raise ValueError if the requirements are incompatible
        with the Index subclass.

        'select' may return an empty object or None if no matches can be
        found.
        """


def select_signature(ss, ksize=None, moltype=None, scaled=0, num=0,
                     containment=False):
    "Check that the given signature matches the specificed requirements."
    # ksize match?
    if ksize and ksize != ss.minhash.ksize:
        return False

    # moltype match?
    if moltype and moltype != ss.minhash.moltype:
        return False

    # containment requires scaled; similarity does not.
    if containment:
        if not scaled:
            raise ValueError("'containment' requires 'scaled' in Index.select'")
        if not ss.minhash.scaled:
            return False

    # 'scaled' and 'num' are incompatible
    if scaled:
        if ss.minhash.num:
            return False
    if num:
        # note, here we check if 'num' is identical; this can be
        # changed later.
        if ss.minhash.scaled or num != ss.minhash.num:
            return False

    return True


class LinearIndex(Index):
    "An Index for a collection of signatures. Can load from a .sig file."
    def __init__(self, _signatures=None, filename=None):
        self._signatures = []
        if _signatures:
            self._signatures = list(_signatures)
        self.filename = filename

    @property
    def location(self):
        return self.filename

    def signatures(self):
        return iter(self._signatures)

    def __bool__(self):
        return bool(self._signatures)

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

    def select(self, **kwargs):
        """Return new LinearIndex containing only signatures that match req's.

        Does not raise ValueError, but may return an empty Index.
        """
        # eliminate things from kwargs with None or zero value
        kw = { k : v for (k, v) in kwargs.items() if v }

        siglist = []
        for ss in self._signatures:
            if select_signature(ss, **kwargs):
                siglist.append(ss)

        return LinearIndex(siglist, self.location)


class LazyLinearIndex(Index):
    "An Index for lazy linear search of another database."
    def __init__(self, db):
        self.db = db

    @property
    def location(self):
        return self.db.location

    def signatures(self):
        for ss in self.db.signatures():
            yield ss

    def __bool__(self):
        try:
            first_sig = next(iter(self.signatures()))
            return True
        except StopIteration:
            return False

    def __len__(self):
        raise NotImplementedError

    def insert(self, node):
        raise NotImplementedError

    def save(self, path):
        raise NotImplementedError

    @classmethod
    def load(cls, path):
        raise NotImplementedError

    def select(self, **kwargs):
        """Return new object yielding only signatures that match req's.

        Does not raise ValueError, but may return an empty Index.
        """
        db = self.db.select(**kwargs)
        return LazyLinearIndex(db)


class ZipFileLinearIndex(Index):
    """\
    A read-only collection of signatures in a zip file.

    Does not support `insert` or `save`.
    """
    is_database = True

    def __init__(self, zf, selection_dict=None,
                 traverse_yield_all=False):
        self.zf = zf
        self.selection_dict = selection_dict
        self.traverse_yield_all = traverse_yield_all

    def __bool__(self):
        # @CTB write test to make sure this doesn't call __len__
        try:
            first_sig = next(iter(self.signatures()))
        except StopIteration:
            return False

        return True

    def __len__(self):
        return len(list(self.signatures()))

    @property
    def location(self):
        return self.zf.filename

    def insert(self, signature):
        raise NotImplementedError

    def save(self, path):
        raise NotImplementedError

    @classmethod
    def load(cls, location, traverse_yield_all=False):
        "Class method to load a zipfile."
        zf = zipfile.ZipFile(location, 'r')
        return cls(zf, traverse_yield_all=traverse_yield_all)

    def signatures(self):
        "Load all signatures in the zip file."
        from .signature import load_signatures
        for zipinfo in self.zf.infolist():
            # should we load this file? if it ends in .sig OR we are forcing:
            if zipinfo.filename.endswith('.sig') or \
               zipinfo.filename.endswith('.sig.gz') or \
               self.traverse_yield_all:
                fp = self.zf.open(zipinfo)

                # now load all the signatures and select on ksize/moltype:
                selection_dict = self.selection_dict

                # note: if 'fp' doesn't contain a valid JSON signature,
                # load_signatures will silently fail & yield nothing.
                for ss in load_signatures(fp):
                    if selection_dict:
                        if select_signature(ss, **self.selection_dict):
                            yield ss
                    else:
                        yield ss

    def select(self, **kwargs):
        "Select signatures in zip file based on ksize/moltype/etc."
        return ZipFileLinearIndex(self.zf,
                                  selection_dict=kwargs,
                                  traverse_yield_all=self.traverse_yield_all)


class QuerySpecific_GatherCounter:
    def __init__(self, query_mh):
        if not query_mh.scaled:
            raise ValueError('gather requires scaled signatures')

        # track query
        self.orig_query_mh = copy.copy(query_mh).flatten()
        self.scaled = query_mh.scaled

        # track matching signatures & their locations
        self.siglist = []
        self.locations = []

        # ...and overlaps with query
        self.counter = Counter()

        # cannot add matches once query has started.
        self.query_started = 0

    def add(self, ss, location=None):
        assert not self.query_started

        i = len(self.siglist)
        self.siglist.append(ss)
        self.locations.append(location)

        # note: scaled will be max of all matches.
        self.downsample(ss.minhash.scaled)

        # upon insertion, count & track overlap with the specific query.
        self.counter[i] = self.orig_query_mh.count_common(ss.minhash, True)

    def downsample(self, scaled):
        if scaled > self.scaled:
            self.scaled = scaled

    def calc_threshold(self, threshold_bp, scaled, query_size):
        # @CTB can be outside this class
        threshold = 0.0
        n_threshold_hashes = 0

        if threshold_bp:
            # if we have a threshold_bp of N, then that amounts to N/scaled
            # hashes:
            n_threshold_hashes = float(threshold_bp) / scaled

            # that then requires the following containment:
            threshold = n_threshold_hashes / query_size

        return threshold, n_threshold_hashes

    def peek(self, cur_query_mh, scaled, threshold_bp=0, **kwargs):
        "Get next potential result."
        self.query_started = 1

        # empty? nothing to search.
        counter = self.counter
        siglist = self.siglist
        if not counter:
            return []

        assert siglist

        self.downsample(scaled)
        scaled = self.scaled
        cur_query_mh = cur_query_mh.downsample(scaled=scaled)

        if not cur_query_mh:             # empty query? quit.
            return []

        assert cur_query_mh.contained_by(self.orig_query_mh,
                                     downsample=True) == 1

        # are we setting a threshold?
        threshold, n_threshold_hashes = self.calc_threshold(threshold_bp,
                                                            scaled,
                                                            len(cur_query_mh))

        # is it too high to ever match? if so, exit.
        if threshold > 1.0:
            return []

        # Find the best match -
        most_common = counter.most_common()
        dataset_id, match_size = most_common[0]

        # below threshold? no match!
        if match_size < n_threshold_hashes:
            return []

        # pull match and location.
        match = siglist[dataset_id]

        # calculate containment
        cont = cur_query_mh.contained_by(match.minhash, downsample=True)

        retval = []
        if cont and cont >= threshold:
            # calculate intersection of this "best match" with query
            # for removal.

            # @CTB: note flatten
            match_mh = match.minhash.downsample(scaled=scaled).flatten()
            intersect_mh = cur_query_mh.intersection(match_mh)
            location = self.locations[dataset_id]
            retval = [IndexSearchResult(cont, match, location), intersect_mh]

        return retval

    def consume(self, intersect_mh):
        "Remove the given hashes."
        siglist = self.siglist
        counter = self.counter

        most_common = counter.most_common()

        # Prepare counter for finding the next match by decrementing
        # all hashes found in the current match in other datasets;
        # remove empty datasets from counter, too.
        for (dataset_id, _) in most_common:
            # @CTB: we may want to downsample remaining_mh here...
            remaining_mh = siglist[dataset_id].minhash
            intersect_count = intersect_mh.count_common(remaining_mh,
                                                        downsample=True)
            counter[dataset_id] -= intersect_count
            if counter[dataset_id] == 0:
                del counter[dataset_id]

    def next(self, query, threshold_bp=0):
        result = self.peek(query.minhash, query.minhash.scaled,
                                   threshold_bp)
        if result:
            (sr, intersect_mh) = result
            self.consume(intersect_mh)
            
            query.minhash.remove_many(intersect_mh.hashes) #  @CTB

            return [sr]
        return []
                


class CounterGatherIndex(Index):
    def __init__(self, query):
        self.counter = QuerySpecific_GatherCounter(query.minhash)

    def insert(self, ss, location=None):
        self.counter.add(ss, location)

    def gather(self, query, threshold_bp=0):
        return self.counter.next(query, threshold_bp)
        
    def signatures(self):
        raise NotImplementedError

    def signatures_with_location(self):
        raise NotImplementedError

    def prefetch(self, *args, **kwargs):
        raise NotImplementedError

    @classmethod
    def load(self, *args):
        raise NotImplementedError

    def save(self, *args):
        raise NotImplementedError

    def find(self, search_fn, *args, **kwargs):
        raise NotImplementedError

    def search(self, query, *args, **kwargs):
        raise NotImplementedError

    def select(self, *args, **kwargs):
        raise NotImplementedError


class MultiIndex(Index):
    """An Index class that wraps other Index classes.

    The MultiIndex constructor takes two arguments: a list of Index
    objects, and a matching list of sources (filenames, etc.)  If the
    source is not None, then it will be used to override the 'filename'
    in the triple that is returned by search and gather.

    One specific use for this is when loading signatures from a directory;
    MultiIndex will properly record which files provided which signatures.
    """
    def __init__(self, index_list, source_list):
        self.index_list = list(index_list)
        self.source_list = list(source_list)
        assert len(index_list) == len(source_list)

    def signatures(self):
        for idx in self.index_list:
            for ss in idx.signatures():
                yield ss

    def signatures_with_location(self):
        for idx, loc in zip(self.index_list, self.source_list):
            for ss in idx.signatures():
                yield ss, loc

    def __len__(self):
        return sum([ len(idx) for idx in self.index_list ])

    def insert(self, *args):
        raise NotImplementedError

    @classmethod
    def load(self, *args):
        raise NotImplementedError

    @classmethod
    def load_from_path(cls, pathname, force=False):
        "Create a MultiIndex from a path (filename or directory)."
        from .sourmash_args import traverse_find_sigs
        if not os.path.exists(pathname): # @CTB change to isdir
            raise ValueError(f"'{pathname}' must be a directory")

        index_list = []
        source_list = []
        for thisfile in traverse_find_sigs([pathname], yield_all_files=force):
            try:
                idx = LinearIndex.load(thisfile)
                index_list.append(idx)
                source_list.append(thisfile)
            except (IOError, sourmash.exceptions.SourmashError):
                if force:
                    continue    # ignore error
                else:
                    raise       # continue past error!

        db = None
        if index_list:
            db = cls(index_list, source_list)
        else:
            raise ValueError(f"no signatures to load under directory '{pathname}'")

        return db

    @classmethod
    def load_from_pathlist(cls, filename):
        "Create a MultiIndex from all files listed in a text file."
        from .sourmash_args import (load_pathlist_from_file,
                                    load_file_as_index)
        idx_list = []
        src_list = []

        file_list = load_pathlist_from_file(filename)
        for fname in file_list:
            idx = load_file_as_index(fname)
            src = fname

            idx_list.append(idx)
            src_list.append(src)

        db = MultiIndex(idx_list, src_list)
        return db

    def save(self, *args):
        raise NotImplementedError

    def select(self, **kwargs):
        "Run 'select' on all indices within this MultiIndex."
        new_idx_list = []
        new_src_list = []
        for idx, src in zip(self.index_list, self.source_list):
            idx = idx.select(**kwargs)
            new_idx_list.append(idx)
            new_src_list.append(src)

        return MultiIndex(new_idx_list, new_src_list)

    def filter(self, filter_fn):
        new_idx_list = []
        new_src_list = []
        for idx, src in zip(self.index_list, self.source_list):
            idx = idx.filter(filter_fn)
            new_idx_list.append(idx)
            new_src_list.append(src)

        return MultiIndex(new_idx_list, new_src_list)

    def search(self, query, **kwargs):
        """Return the match with the best Jaccard similarity in the Index.

        Note: this overrides the location of the match if needed.
        """
        # do the actual search:
        matches = []
        for idx, src in zip(self.index_list, self.source_list):
            for (score, ss, filename) in idx.search(query, **kwargs):
                best_src = src or filename # override if src provided
                matches.append(IndexSearchResult(score, ss, best_src))
                
        # sort!
        matches.sort(key=lambda x: -x.score)
        return matches

    def prefetch(self, query, threshold_bp, **kwargs):
        "Return all matches with specified overlap."
        # actually do search!
        results = []
        for idx, src in zip(self.index_list, self.source_list):
            if not idx:
                continue

            for (score, ss, filename) in idx.prefetch(query, threshold_bp,
                                                      **kwargs):
                best_src = src or filename # override if src provided
                yield IndexSearchResult(score, ss, best_src)
            
        return results

    # note: 'gather' is inherited from Index base class, and uses prefetch.
