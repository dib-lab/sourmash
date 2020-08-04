#! /usr/bin/env python
"""
Walk through a file of long contigs and analyze each one for contamination
against an SBT database.

See http://ivory.idyll.org/blog/2018-detecting-contamination-in-long-read-assemblies.html as well.

Author: C. Titus Brown, titus@idyll.org

Requires sourmash 3.x and Python 3.6
"""
import sys
import argparse
import screed
import sourmash
import csv


def main():
    p = argparse.ArgumentParser()
    p.add_argument('input_seqs')
    p.add_argument('sbt_database')
    p.add_argument('--output-nomatch', type=argparse.FileType('wt'))
    p.add_argument('--output-match', type=argparse.FileType('wt'))
    p.add_argument('--csv', type=argparse.FileType('wt'))
    p.add_argument('-t', '--use-tabs', action='store_true')
    p.add_argument('--threshold-bp', default=10000, type=float,
                   help='est number of bp to require before reporting match')
    args = p.parse_args()

    tree = sourmash.load_sbt_index(args.sbt_database)
    print(f'found SBT database {args.sbt_database}')

    leaf = next(iter(tree.leaves()))
    mh = leaf.data.minhash.copy_and_clear()

    print(f'using ksize={mh.ksize}, scaled={mh.scaled}')

    print(f'loading sequences from {args.input_seqs}')
    if args.output_match:
        print(f'saving match sequences to {args.output_match.name}')
    if args.output_nomatch:
        print(f'saving nomatch sequences to {args.output_nomatch.name}')
    if args.csv:
        print(f'outputting CSV summary to {args.csv.name}')

    found_list = []
    total_seqs = 0
    total_bp = 0
    match_seqs = 0
    match_bp = 0
    nomatch_seqs = 0
    nomatch_bp = 0
    
    matches = 0
    for record in screed.open(args.input_seqs):
        total_seqs += 1
        total_bp += len(record.sequence)
        found = False

        query_mh = mh.copy_and_clear()
        query_mh.add_sequence(record.sequence, force=True)
        query = sourmash.SourmashSignature(query_mh)

        # too small a sequence/not enough hashes? notify
        if not query_mh.get_mins():
            print(f'note: skipping {query.name()[:20]}, no hashes in sketch')
            continue

        for result in tree.gather(query, threshold_bp=query_mh.scaled * len(query_mh)):
            (containment, match, name) = result

            in_common = match.minhash.count_common(query_mh)
            contamination_min = in_common * query_mh.scaled
            if contamination_min < args.threshold_bp:
                break

            found = True
            matches += 1
            found_list.append((record.name, match.name(), contamination_min))
            break

        if not found:
            found_list.append((record.name, '', 0.0))

        if found:
            if args.output_match:
                args.output_match.write(f'>{record.name}\n{record.sequence}\n')
            match_seqs += 1
            match_bp += len(record.sequence)
        else: # not found
            if args.output_nomatch:
                args.output_nomatch.write(f'>{record.name}\n{record.sequence}\n')
            nomatch_seqs += 1
            nomatch_bp += len(record.sequence)

        print(f'PROGRESS: searched {total_seqs}, found {matches} ({int(matches/total_seqs*100)}%)')
        if found:
            print(f'found: est {int(contamination_min/1e3)} kb matches to {match.name()}')
        sys.stdout.flush()

    print(f'searched {total_seqs} seqs with {int(total_bp/1e3)}kb, found {matches} in database')
    print(f'MATCHED {match_seqs} sequences with {int(match_bp/1e3)}kb.')
    print(f'NOMATCH {nomatch_seqs} sequences with {int(nomatch_bp/1e3)}kb.')

    print('')

    if args.csv:
        delimiter = ','
        if args.use_tabs:
            delimiter = '\t'
        w = csv.DictWriter(args.csv,
                           fieldnames=['query', 'match', 'min_bp_found'],
                           delimiter=delimiter)
        w.writeheader()
        for (query, match, min_bp) in found_list:
            w.writerow(dict(query=query, match=match, min_bp_found=min_bp))

if __name__ == '__main__':
    main()
