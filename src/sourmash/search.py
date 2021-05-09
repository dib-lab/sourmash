"""
Code for searching collections of signatures.
"""
from collections import namedtuple
import sys
import os
from enum import Enum

from .logging import notify, error
from .signature import SourmashSignature
from .minhash import _get_max_hash_for_scaled


class SearchType(Enum):
    JACCARD = 1
    CONTAINMENT = 2
    MAX_CONTAINMENT = 3


def make_jaccard_search_query(*,
                              do_containment=False,
                              do_max_containment=False,
                              best_only=False,
                              threshold=None):
    """\
    Make a "flat" search object for Jaccard search & containment.
    """
    if do_containment and do_max_containment:
        raise TypeError("'do_containment' and 'do_max_containment' cannot both be True")

    # configure search - containment? ignore abundance? best only?
    search_cls = JaccardSearch
    if best_only:
        search_cls = JaccardSearchBestOnly

    if do_containment:
        search_obj = search_cls(SearchType.CONTAINMENT, threshold)
    elif do_max_containment:
        search_obj = search_cls(SearchType.MAX_CONTAINMENT, threshold)
    else:
        search_obj = search_cls(SearchType.JACCARD, threshold)

    return search_obj


def make_gather_query(query_mh, threshold_bp, *, best_only=True):
    "Make a search object for gather."
    if not query_mh:
        raise ValueError("query is empty!?")

    scaled = query_mh.scaled
    if not scaled:
        raise TypeError("query signature must be calculated with scaled")

    # are we setting a threshold?
    threshold = 0
    if threshold_bp:
        if threshold_bp < 0:
            raise TypeError("threshold_bp must be non-negative")

        # if we have a threshold_bp of N, then that amounts to N/scaled
        # hashes:
        n_threshold_hashes = threshold_bp / scaled

        # that then requires the following containment:
        threshold = n_threshold_hashes / len(query_mh)

        # is it too high to ever match? if so, exit.
        if threshold > 1.0:
            raise ValueError("requested threshold_bp is unattainable with this query")

    if best_only:
        search_obj = JaccardSearchBestOnly(SearchType.CONTAINMENT,
                                           threshold=threshold)
    else:
        search_obj = JaccardSearch(SearchType.CONTAINMENT,
                                   threshold=threshold)

    return search_obj


class JaccardSearch:
    """
    A class used by Index classes for searching/gathering.
    """
    def __init__(self, search_type, threshold=None):
        "Constructor. Takes type of search, and optional threshold."
        score_fn = None
        require_scaled = False

        if search_type == SearchType.JACCARD:
            score_fn = self.score_jaccard
        elif search_type == SearchType.CONTAINMENT:
            score_fn = self.score_containment
            require_scaled = True
        elif search_type == SearchType.MAX_CONTAINMENT:
            score_fn = self.score_max_containment
            require_scaled = True
        self.score_fn = score_fn
        self.require_scaled = require_scaled

        if threshold is None:
            threshold = 0
        self.threshold = float(threshold)

    def check_is_compatible(self, sig):
        """
        Is this query compatible with this type of search? Raise TypeError
        if not.
        """
        if self.require_scaled:
            if not sig.minhash.scaled:
                raise TypeError("this search requires a scaled signature")

        if sig.minhash.track_abundance:
            raise TypeError("this search cannot be done with an abund signature")

    def passes(self, score):
        """Return True if this score meets or exceeds the threshold.

        Note: this can be used whenever a score or estimate is available
        (e.g. internal nodes on an SBT). `collect(...)`, below, decides
        whether a particular signature should be collected, and/or can
        update the threshold (used for BestOnly behavior).
        """
        if score and score >= self.threshold:
            return True
        return False

    def collect(self, score, match_sig):
        "Return True if this match should be collected."
        return True

    def score_jaccard(self, query_size, shared_size, subject_size, total_size):
        "Calculate Jaccard similarity."
        return shared_size / total_size

    def score_containment(self, query_size, shared_size, subject_size,
                          total_size):
        "Calculate Jaccard containment."
        if query_size == 0:
            return 0
        return shared_size / query_size

    def score_max_containment(self, query_size, shared_size, subject_size,
                              total_size):
        "Calculate Jaccard max containment."
        min_denom = min(query_size, subject_size)
        if min_denom == 0:
            return 0
        return shared_size / min_denom


class JaccardSearchBestOnly(JaccardSearch):
    "A subclass of JaccardSearch that implements best-only."
    def collect(self, score, match):
        "Raise the threshold to the best match found so far."
        self.threshold = max(self.threshold, score)
        return True


# generic SearchResult tuple.
SearchResult = namedtuple('SearchResult',
                          'similarity, match, md5, filename, name, query, query_filename, query_name, query_md5')


def format_bp(bp):
    "Pretty-print bp information."
    bp = float(bp)
    if bp < 500:
        return '{:.0f} bp '.format(bp)
    elif bp <= 500e3:
        return '{:.1f} kbp'.format(round(bp / 1e3, 1))
    elif bp < 500e6:
        return '{:.1f} Mbp'.format(round(bp / 1e6, 1))
    elif bp < 500e9:
        return '{:.1f} Gbp'.format(round(bp / 1e9, 1))
    return '???'


def search_databases_with_flat_query(query, databases, **kwargs):
    results = []
    found_md5 = set()

    for db in databases:
        search_iter = db.search(query, **kwargs)
        for (score, match, filename) in search_iter:
            md5 = match.md5sum()
            if md5 not in found_md5:
                results.append((score, match, filename))
                found_md5.add(md5)

    # sort results on similarity (reverse)
    results.sort(key=lambda x: -x[0])

    x = []
    for (score, match, filename) in results:
        x.append(SearchResult(similarity=score,
                              match=match,
                              md5=match.md5sum(),
                              filename=filename,
                              name=match.name,
                              query=query,
                              query_filename=query.filename,
                              query_name=query.name,
                              query_md5=query.md5sum()[:8]
        ))
    return x


def search_databases_with_abund_query(query, databases, **kwargs):
    results = []
    found_md5 = set()

    for db in databases:
        search_iter = db.search_abund(query, **kwargs)
        for (score, match, filename) in search_iter:
            md5 = match.md5sum()
            if md5 not in found_md5:
                results.append((score, match, filename))
                found_md5.add(md5)

    # sort results on similarity (reverse)
    results.sort(key=lambda x: -x[0])

    x = []
    for (score, match, filename) in results:
        x.append(SearchResult(similarity=score,
                              match=match,
                              md5=match.md5sum(),
                              filename=filename,
                              name=match.name,
                              query=query,
                              query_filename=query.filename,
                              query_name=query.name,
                              query_md5=query.md5sum()[:8]
        ))
    return x

###
### gather code
###

GatherResult = namedtuple('GatherResult',
                          'intersect_bp, f_orig_query, f_match, f_unique_to_query, f_unique_weighted, average_abund, median_abund, std_abund, filename, name, md5, match, f_match_orig, unique_intersect_bp, gather_result_rank, remaining_bp')


# build a new query object, subtracting found mins and downsampling
def _subtract_and_downsample(to_remove, old_query, scaled=None):
    mh = old_query.minhash.mutable()
    mh = mh.downsample(scaled=scaled)
    mh.remove_many(to_remove)

    return SourmashSignature(mh)


def _find_best(counters, query, threshold_bp):
    """
    Search for the best containment, return precisely one match.
    """
    results = []

    best_result = None
    best_intersect_mh = None

    # find the best score across multiple counters, without consuming
    for counter in counters:
        result = counter.peek(query.minhash, threshold_bp)
        if result:
            (sr, intersect_mh) = result

            if best_result is None or sr.score > best_result.score:
                best_result = sr
                best_intersect_mh = intersect_mh

    if best_result:
        # remove the best result from each counter
        for counter in counters:
            counter.consume(best_intersect_mh)

        # and done!
        return best_result
    return None


def _filter_max_hash(values, max_hash):
    results = set()
    for v in values:
        if v < max_hash:
            results.add(v)
    return results


def gather_databases(query, counters, threshold_bp, ignore_abundance):
    """
    Iteratively find the best containment of `query` in all the `counters`,
    until we find fewer than `threshold_bp` (estimated) bp in common.
    """
    # track original query information for later usage.
    track_abundance = query.minhash.track_abundance and not ignore_abundance
    orig_query_mh = query.minhash
    orig_query_hashes = set(orig_query_mh.hashes)

    # do we pay attention to abundances?
    orig_query_abunds = { k: 1 for k in orig_query_hashes }
    if track_abundance:
        import numpy as np
        orig_query_abunds = orig_query_mh.hashes

    query.minhash = query.minhash.flatten()

    cmp_scaled = query.minhash.scaled    # initialize with resolution of query
    result_n = 0
    while query.minhash:
        # find the best match!
        best_result = _find_best(counters, query, threshold_bp)

        if not best_result:          # no matches at all for this cutoff!
            notify(f'found less than {format_bp(threshold_bp)} in common. => exiting')
            break

        best_match = best_result.signature
        filename = best_result.location

        # subtract found hashes from search hashes, construct new search
        query_hashes = set(query.minhash.hashes)
        found_hashes = set(best_match.minhash.hashes)

        # Is the best match computed with scaled? Die if not.
        match_scaled = best_match.minhash.scaled
        assert match_scaled

        # pick the highest scaled / lowest resolution
        cmp_scaled = max(cmp_scaled, match_scaled)

        # eliminate hashes under this new resolution.
        # (CTB note: this means that if a high scaled/low res signature is
        # found early on, resolution will be low from then on.)
        new_max_hash = _get_max_hash_for_scaled(cmp_scaled)
        query_hashes = _filter_max_hash(query_hashes, new_max_hash)
        found_hashes = _filter_max_hash(found_hashes, new_max_hash)
        orig_query_hashes = _filter_max_hash(orig_query_hashes, new_max_hash)
        sum_abunds = sum(( orig_query_abunds[k] for k in orig_query_hashes))

        # calculate intersection with query hashes:
        intersect_hashes = query_hashes.intersection(found_hashes)
        unique_intersect_bp = cmp_scaled * len(intersect_hashes)
        intersect_orig_hashes = orig_query_hashes.intersection(found_hashes)
        intersect_bp = cmp_scaled * len(intersect_orig_hashes)

        # calculate fractions wrt first denominator - genome size
        assert intersect_hashes.issubset(found_hashes)
        f_match = len(intersect_hashes) / len(found_hashes)
        f_orig_query = len(intersect_orig_hashes) / len(orig_query_hashes)

        # calculate fractions wrt second denominator - metagenome size
        assert intersect_hashes.issubset(orig_query_hashes)
        f_unique_to_query = len(intersect_hashes) / len(orig_query_hashes)

        # calculate fraction of subject match with orig query
        f_match_orig = best_match.minhash.contained_by(orig_query_mh,
                                                       downsample=True)

        # calculate scores weighted by abundances
        f_unique_weighted = sum((orig_query_abunds[k] for k in intersect_hashes))
        f_unique_weighted /= sum_abunds

        # calculate stats on abundances, if desired.
        average_abund, median_abund, std_abund = None, None, None
        if track_abundance:
            intersect_abunds = (orig_query_abunds[k] for k in intersect_hashes)
            intersect_abunds = list(intersect_abunds)

            average_abund = np.mean(intersect_abunds)
            median_abund = np.median(intersect_abunds)
            std_abund = np.std(intersect_abunds)

        # construct a new query, subtracting hashes found in previous one.
        query = _subtract_and_downsample(found_hashes, query, cmp_scaled)
        remaining_bp = cmp_scaled * len(query.minhash.hashes)

        # compute weighted_missed:
        query_hashes -= set(found_hashes)
        weighted_missed = sum((orig_query_abunds[k] for k in query_hashes)) \
             / sum_abunds

        # build a result namedtuple
        result = GatherResult(intersect_bp=intersect_bp,
                              unique_intersect_bp=unique_intersect_bp,
                              f_orig_query=f_orig_query,
                              f_match=f_match,
                              f_match_orig=f_match_orig,
                              f_unique_to_query=f_unique_to_query,
                              f_unique_weighted=f_unique_weighted,
                              average_abund=average_abund,
                              median_abund=median_abund,
                              std_abund=std_abund,
                              filename=filename,
                              md5=best_match.md5sum(),
                              name=str(best_match),
                              match=best_match,
                              gather_result_rank=result_n,
                              remaining_bp=remaining_bp)
        result_n += 1

        yield result, weighted_missed, new_max_hash, query


###
### prefetch code
###

PrefetchResult = namedtuple('PrefetchResult',
                            'intersect_bp, jaccard, max_containment, f_query_match, f_match_query, match, match_filename, match_name, match_md5, match_bp, query, query_filename, query_name, query_md5, query_bp')


def prefetch_database(query, database, threshold_bp):
    """
    Find all matches to `query_mh` >= `threshold_bp` in `database`.
    """
    query_mh = query.minhash
    scaled = query_mh.scaled
    assert scaled

    # for testing/double-checking purposes, calculate expected threshold -
    threshold = threshold_bp / scaled

    # iterate over all signatures in database, find matches

    for result in database.prefetch(query, threshold_bp):
        # base intersections on downsampled minhashes
        match = result.signature
        db_mh = match.minhash.flatten().downsample(scaled=scaled)

        # calculate db match intersection with query hashes:
        intersect_mh = query_mh.intersection(db_mh)
        assert len(intersect_mh) >= threshold

        f_query_match = db_mh.contained_by(query_mh)
        f_match_query = query_mh.contained_by(db_mh)
        max_containment = max(f_query_match, f_match_query)

        # build a result namedtuple
        result = PrefetchResult(
            intersect_bp=len(intersect_mh) * scaled,
            query_bp = len(query_mh) * scaled,
            match_bp = len(db_mh) * scaled,
            jaccard=db_mh.jaccard(query_mh),
            max_containment=max_containment,
            f_query_match=f_query_match,
            f_match_query=f_match_query,
            match=match,
            match_filename=match.filename,
            match_name=match.name,
            match_md5=match.md5sum()[:8],
            query=query,
            query_filename=query.filename,
            query_name=query.name,
            query_md5=query.md5sum()[:8]
        )

        yield result
