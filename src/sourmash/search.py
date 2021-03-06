from collections import namedtuple
import sys

from .logging import notify, error
from .signature import SourmashSignature
from .minhash import _get_max_hash_for_scaled


# generic SearchResult.
SearchResult = namedtuple('SearchResult',
                          'similarity, match, md5, filename, name')


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


def search_databases(query, databases, threshold, do_containment, best_only,
                     ignore_abundance, unload_data=False):
    results = []
    found_md5 = set()
    for (obj, filename, filetype) in databases:
        search_iter = obj.search(query, threshold=threshold,
                                 do_containment=do_containment,
                                 ignore_abundance=ignore_abundance,
                                 best_only=best_only,
                                 unload_data=unload_data)
        for (similarity, match, filename) in search_iter:
            md5 = match.md5sum()
            if md5 not in found_md5:
                results.append((similarity, match, filename))
                found_md5.add(md5)

    # sort results on similarity (reverse)
    results.sort(key=lambda x: -x[0])

    x = []
    for (similarity, match, filename) in results:
        x.append(SearchResult(similarity=similarity,
                              match=match,
                              md5=match.md5sum(),
                              filename=filename,
                              name=match.name))
    return x

###
### gather code
###

GatherResult = namedtuple('GatherResult',
                          'intersect_bp, f_orig_query, f_match, f_unique_to_query, f_unique_weighted, average_abund, median_abund, std_abund, filename, name, md5, match, f_match_orig, unique_intersect_bp, gather_result_rank, remaining_bp')


# build a new query object, subtracting found mins and downsampling
def _subtract_and_downsample(to_remove, old_query, scaled=None):
    mh = old_query.minhash
    mh = mh.downsample(scaled=scaled)
    mh.remove_many(to_remove)

    return SourmashSignature(mh)


def _find_best(dblist, query, threshold_bp):
    """
    Search for the best containment, return precisely one match.
    """

    best_cont = 0.0
    best_match = None
    best_filename = None

    # quantize threshold_bp to be an integer multiple of scaled
    query_scaled = query.minhash.scaled
    threshold_bp = int(threshold_bp / query_scaled) * query_scaled

    # search across all databases
    for (obj, filename, filetype) in dblist:
        for cont, match, fname in obj.counter_gather(query, threshold_bp=threshold_bp):
            assert cont                   # all matches should be nonzero.

            # note, break ties based on name, to ensure consistent order.
            if (cont == best_cont and str(match) < str(best_match)) or \
               cont > best_cont:
                # update best match.
                best_cont = cont
                best_match = match

                # some objects may not have associated filename (e.g. SBTs)
                best_filename = fname or filename

    if not best_match:
        return None, None, None

    return best_cont, best_match, best_filename


def _filter_max_hash(values, max_hash):
    results = set()
    for v in values:
        if v < max_hash:
            results.add(v)
    return results


def gather_databases(query, databases, threshold_bp, ignore_abundance):
    """
    Iteratively find the best containment of `query` in all the `databases`,
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

    cmp_scaled = query.minhash.scaled    # initialize with resolution of query
    result_n = 0
    while query.minhash:
        # find the best match!
        best_cont, best_match, filename = _find_best(databases, query,
                                                     threshold_bp)
        if not best_match:          # no matches at all for this cutoff!
            notify('found less than {} in common. => exiting',
                   format_bp(threshold_bp))
            break

        # subtract found hashes from search hashes, construct new search
        query_hashes = set(query.minhash.hashes)
        found_hashes = set(best_match.minhash.hashes)

        # Is the best match computed with scaled? Die if not.
        match_scaled = best_match.minhash.scaled
        if not match_scaled:
            error('Best match in gather is not scaled.')
            error('Please prepare gather databases with --scaled')
            raise Exception

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
