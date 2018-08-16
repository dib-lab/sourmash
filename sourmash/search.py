from __future__ import division
from collections import namedtuple
import sys

from .logging import notify, error
from .signature import SourmashSignature
from .sbtmh import (search_minhashes,
                    search_minhashes_containment,
                    select_signature)
from .sbtmh import SearchMinHashesFindBest, SearchMinHashesFindBestIgnoreMaxHash
from ._minhash import get_max_hash_for_scaled


# generic SearchResult across individual signatures + SBTs.
SearchResult = namedtuple('SearchResult',
                          'similarity, match_sig, md5, filename, name')


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


def search_databases(query, databases, threshold, do_containment, best_only):
    # set up the search & score function(s) - similarity vs containment
    search_fn = search_minhashes
    query_match = lambda x: query.similarity(x, downsample=True)
    if do_containment:
        search_fn = search_minhashes_containment
        query_match = lambda x: query.contained_by(x, downsample=True)

    results = []
    found_md5 = set()
    for (sbt_or_siglist, filename, is_sbt) in databases:
        if is_sbt:
            if best_only:            # this needs to be reset for each SBT
                search_fn = SearchMinHashesFindBest().search

            tree = sbt_or_siglist
            for leaf in tree.find(search_fn, query, threshold):
                leafdata = select_signature(leaf, query)
                similarity = query_match(leafdata)

                # tree search should always/only return matches above threshold
                assert similarity >= threshold

                if leafdata.md5sum() not in found_md5:
                    sr = SearchResult(similarity=similarity,
                                      match_sig=leafdata,
                                      md5=leafdata.md5sum(),
                                      filename=filename,
                                      name=leafdata.name())
                    found_md5.add(sr.md5)
                    results.append(sr)

        else: # list of signatures
            for ss in sbt_or_siglist:
                similarity = query_match(ss)
                if similarity >= threshold and \
                       ss.md5sum() not in found_md5:
                    sr = SearchResult(similarity=similarity,
                                      match_sig=ss,
                                      md5=ss.md5sum(),
                                      filename=filename,
                                      name=ss.name())
                    found_md5.add(sr.md5)
                    results.append(sr)


    # sort results on similarity (reverse)
    results.sort(key=lambda x: -x.similarity)

    return results


GatherResult = namedtuple('GatherResult',
                          'intersect_bp, f_orig_query, f_match, f_unique_to_query, f_unique_weighted, average_abund, median_abund, std_abund, filename, name, md5, leaf')

def gather_databases(query, databases, threshold_bp, ignore_abundance):
    orig_query = query
    orig_mins = orig_query.minhash.get_hashes()
    orig_abunds = { k: 1 for k in orig_mins }

    # do we pay attention to abundances?
    if orig_query.minhash.track_abundance and not ignore_abundance:
        import numpy as np
        orig_abunds = orig_query.minhash.get_mins(with_abundance=True)

    # calculate the band size/resolution R for the genome
    R_metagenome = orig_query.minhash.scaled

    # define a function to do a 'best' search and get only top match.
    def find_best(dblist, query):
        # CTB: could optimize by sharing scores across searches, i.e.
        # a good early score truncates later searches.

        results = []
        for (sbt_or_siglist, filename, is_sbt) in dblist:
            # search a tree
            if is_sbt:
                tree = sbt_or_siglist
                search_fn = SearchMinHashesFindBestIgnoreMaxHash().search

                for leaf in tree.find(search_fn, query, 0.0):
                    leafdata = select_signature(leaf, query)
                    leaf_e = leafdata.minhash
                    similarity = query.minhash.similarity_ignore_maxhash(leaf_e)
                    if similarity > 0.0:
                        results.append((similarity, leafdata))

            # search a signature
            else:
                for ss in sbt_or_siglist:
                    similarity = query.minhash.similarity_ignore_maxhash(ss.minhash)
                    if similarity > 0.0:
                        results.append((similarity, ss))

        if not results:
            return None, None, None

        # take the best result
        results.sort(key=lambda x: (-x[0], x[1].name()))   # reverse sort on similarity, and then on name
        best_similarity, best_leaf = results[0]
        return best_similarity, best_leaf, filename


    # define a function to build new signature object from set of mins
    def build_new_signature(mins, template_sig, scaled=None):
        e = template_sig.minhash.copy_and_clear()
        e.add_many(mins)
        if scaled:
            e = e.downsample_scaled(scaled)
        return SourmashSignature(e)

    # construct a new query that doesn't have the max_hash attribute set.
    new_mins = query.minhash.get_hashes()
    query = build_new_signature(new_mins, orig_query)

    R_comparison = 0
    while 1:
        best_similarity, best_leaf, filename = find_best(databases, query)
        if not best_leaf:          # no matches at all!
            break

        # subtract found hashes from search hashes, construct new search
        query_mins = set(query.minhash.get_hashes())
        found_mins = best_leaf.minhash.get_hashes()

        # figure out what the resolution of the banding on the subject is
        if not best_leaf.minhash.max_hash:
            error('Best hash match in sbt_gather has no max_hash')
            error('Please prepare database of sequences with --scaled')
            sys.exit(-1)

        R_genome = best_leaf.minhash.scaled

        # pick the highest R / lowest resolution
        R_comparison = max(R_comparison, R_metagenome, R_genome)

        # eliminate mins under this new resolution.
        # (CTB note: this means that if a high scaled/low res signature is
        # found early on, resolution will be low from then on.)
        new_max_hash = get_max_hash_for_scaled(R_comparison)
        query_mins = set([ i for i in query_mins if i < new_max_hash ])
        found_mins = set([ i for i in found_mins if i < new_max_hash ])
        orig_mins = set([ i for i in orig_mins if i < new_max_hash ])
        sum_abunds = sum([ v for (k,v) in orig_abunds.items() if k < new_max_hash ])

        # calculate intersection:
        intersect_mins = query_mins.intersection(found_mins)
        intersect_orig_mins = orig_mins.intersection(found_mins)
        intersect_bp = R_comparison * len(intersect_orig_mins)

        if intersect_bp < threshold_bp:   # hard cutoff for now
            notify('found less than {} in common. => exiting',
                   format_bp(intersect_bp))
            break

        # calculate fractions wrt first denominator - genome size
        genome_n_mins = len(found_mins)
        f_match = len(intersect_mins) / float(genome_n_mins)
        f_orig_query = len(intersect_orig_mins) / float(len(orig_mins))

        # calculate fractions wrt second denominator - metagenome size
        orig_mh = orig_query.minhash.downsample_scaled(R_comparison)
        query_n_mins = len(orig_mh)
        f_unique_to_query = len(intersect_mins) / float(query_n_mins)

        # calculate scores weighted by abundances
        f_unique_weighted = sum((orig_abunds[k] for k in intersect_mins)) \
               / sum_abunds

        intersect_abunds = list(sorted(orig_abunds[k] for k in intersect_mins))
        average_abund, median_abund, std_abund = 0, 0, 0
        if orig_query.minhash.track_abundance and not ignore_abundance:
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
                              md5=best_leaf.md5sum(),
                              name=best_leaf.name(),
                              leaf=best_leaf)

        # construct a new query, minus the previous one.
        query_mins -= set(found_mins)
        query = build_new_signature(query_mins, orig_query, R_comparison)

        weighted_missed = sum((orig_abunds[k] for k in query_mins)) \
             / sum_abunds

        yield result, weighted_missed, new_max_hash, query
