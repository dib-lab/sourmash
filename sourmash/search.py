from __future__ import division
from collections import namedtuple
import sys

from .logging import notify, error
from .signature import SourmashSignature
from ._minhash import get_max_hash_for_scaled


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
                              name=match.name()))
    return x

###
### gather code
###

GatherResult = namedtuple('GatherResult',
                          'intersect_bp, f_orig_query, f_match, f_unique_to_query, f_unique_weighted, average_abund, median_abund, std_abund, filename, name, md5, match')


# build a new query object, subtracting found mins and downsampling
def _subtract_and_downsample(to_remove, old_query, scaled=None):
    mh = old_query.minhash
    mh = mh.downsample_scaled(scaled)
    mh.remove_many(to_remove)

    return SourmashSignature(mh)


def _find_best(dblist, query, threshold_bp):
    """
    Search for the best containment, return precisely one match.
    """

    best_cont = 0.0
    best_match = None
    best_filename = None

    # search across all databases
    # @CTB filetype no longer needed here
    for (obj, filename, filetype) in dblist:
        for cont, match, fname in obj.gather(query, threshold_bp=threshold_bp):
            assert cont                   # all matches should be nonzero.

            # note, break ties based on name, to ensure consistent order.
            if (cont == best_cont and match.name() < best_match.name()) or \
               cont > best_cont:
                # update best match.
                best_cont = cont
                best_match = match

                # some objects may not have associated filename (e.g. SBTs)
                best_filename = fname or filename

    if not best_match:
        return None, None, None

    return best_cont, best_match, best_filename


def gather_databases(query, databases, threshold_bp, ignore_abundance):
    """
    Iteratively find the best containment of `query` in all the `databases`,
    until we find fewer than `threshold_bp` (estimated) bp in common.
    """
    # track original query information for later usage.
    track_abundance = query.minhash.track_abundance and not ignore_abundance
    orig_mh = query.minhash
    orig_mins = orig_mh.get_hashes()
    orig_abunds = { k: 1 for k in orig_mins }

    # do we pay attention to abundances?
    if track_abundance:
        import numpy as np
        orig_abunds = orig_mh.get_mins(with_abundance=True)

    cmp_scaled = query.minhash.scaled    # initialize with resolution of query
    while query.minhash:
        best_cont, best_match, filename = _find_best(databases, query,
                                                     threshold_bp)
        if not best_match:          # no matches at all!
            break

        # subtract found hashes from search hashes, construct new search
        query_mins = set(query.minhash.get_hashes())
        found_mins = best_match.minhash.get_hashes()

        # Is the best match computed with scaled? Die if not.
        match_scaled = best_match.minhash.scaled
        if not match_scaled:
            error('Best match in gather is not scaled.')
            error('Please prepare gather databases with --scaled')
            raise Exception

        # pick the highest scaled / lowest resolution
        cmp_scaled = max(cmp_scaled, match_scaled)

        # eliminate mins under this new resolution.
        # (CTB note: this means that if a high scaled/low res signature is
        # found early on, resolution will be low from then on.)
        new_max_hash = get_max_hash_for_scaled(cmp_scaled)
        query_mins = set([ i for i in query_mins if i < new_max_hash ])
        found_mins = set([ i for i in found_mins if i < new_max_hash ])
        orig_mins = set([ i for i in orig_mins if i < new_max_hash ])
        sum_abunds = sum([ v for (k,v) in orig_abunds.items() if k < new_max_hash ])

        # calculate intersection:
        intersect_mins = query_mins.intersection(found_mins)
        intersect_orig_mins = orig_mins.intersection(found_mins)
        intersect_bp = cmp_scaled * len(intersect_orig_mins)

        if intersect_bp < threshold_bp:   # hard cutoff for now
            notify('found less than {} in common. => exiting',
                   format_bp(threshold_bp))
            break

        # calculate fractions wrt first denominator - genome size
        genome_n_mins = len(found_mins)
        f_match = len(intersect_mins) / float(genome_n_mins)
        f_orig_query = len(intersect_orig_mins) / float(len(orig_mins))

        # calculate fractions wrt second denominator - metagenome size
        orig_mh = orig_mh.downsample_scaled(cmp_scaled)
        query_n_mins = len(orig_mh)
        f_unique_to_query = len(intersect_mins) / float(query_n_mins)

        # calculate scores weighted by abundances
        f_unique_weighted = sum((orig_abunds[k] for k in intersect_mins)) \
               / sum_abunds

        # calculate stats on abundances, if desired.
        average_abund, median_abund, std_abund = 0, 0, 0
        if track_abundance:
            intersect_abunds = list((orig_abunds[k] for k in intersect_mins))
            average_abund = np.mean(intersect_abunds)
            median_abund = np.median(intersect_abunds)
            std_abund = np.std(intersect_abunds)

        # build a result namedtuple
        result = GatherResult(intersect_bp=intersect_bp,
                              f_orig_query=f_orig_query,
                              f_match=f_match,
                              f_unique_to_query=f_unique_to_query,
                              f_unique_weighted=f_unique_weighted,
                              average_abund=average_abund,
                              median_abund=median_abund,
                              std_abund=std_abund,
                              filename=filename,
                              md5=best_match.md5sum(),
                              name=best_match.name(),
                              match=best_match)

        # construct a new query, subtracting hashes found in previous one.
        query = _subtract_and_downsample(found_mins, query, cmp_scaled)

        # compute weighted_missed:
        query_mins -= set(found_mins)
        weighted_missed = sum((orig_abunds[k] for k in query_mins)) \
             / sum_abunds

        yield result, weighted_missed, new_max_hash, query
