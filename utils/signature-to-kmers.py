#! /usr/bin/env python3
"""
Given a signature file and a collection of sequences, output all of the
k-mers and sequences that match a hashval in the signature file.

NOTE: for now, only works for DNA.
"""
import sys
import argparse
import sourmash
from sourmash import MinHash
from sourmash import sourmash_args
from sourmash.minhash import hash_murmur
import screed
import csv
from sourmash.logging import notify, error


NOTIFY_EVERY_BP=1e5
CHUNKSIZE=1e5


def iterate_sequence_chunks(record, ksize):
    global CHUNKSIZE
    CHUNKSIZE = int(CHUNKSIZE)

    for start in range(0, len(record.sequence), CHUNKSIZE):
        yield record.name, record.sequence[start:start + CHUNKSIZE + ksize]


def get_kmers_for_hashvals(sequence, hashvals, ksize):
    "Return k-mers from 'sequence' that yield hashes in 'hashvals'."
    # uppercase!
    sequence = sequence.upper()

    for start in range(0, len(sequence) - ksize + 1):
        kmer = sequence[start:start + ksize]
        kmer_rc = screed.rc(kmer)
        if kmer > kmer_rc:                # choose fwd or rc
            kmer = kmer_rc

        # NOTE: we do not avoid non-ACGT characters, because those k-mers,
        # when hashed, shouldn't match anything that sourmash outputs.
        hashval = hash_murmur(kmer)
        if hashval in hashvals:
            yield kmer, hashval


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query') 					# signature file
    p.add_argument('seqfiles', nargs='+')		# sequence files from which to look for matches
    p.add_argument('--output-sequences', type=str, default=None,
                   help='save matching sequences to this file.')
    p.add_argument('--output-kmers', type=str, default=None,
                   help='save matching kmers to this file.')
    p.add_argument('-k', '--ksize', default=31)
    args = p.parse_args()

    # set up the outputs.
    seqout_fp = None
    if args.output_sequences:
        seqout_fp = open(args.output_sequences, 'wt')

    kmerout_fp = None
    if args.output_kmers:
        kmerout_fp = open(args.output_kmers, 'wt')
        kmerout_w = csv.writer(kmerout_fp)
        kmerout_w.writerow(['kmer', 'hashval'])

    if not (seqout_fp or kmerout_fp):
        error("No output options given!")
        return(-1)

    # first, load the signature and extract the hashvals
    sigobj = sourmash.load_one_signature(args.query, args.ksize)
    assert sigobj.minhash.moltype == 'DNA'
    query_hashvals = set(sigobj.minhash.hashes)
    query_ksize = sigobj.minhash.ksize
    new_mh = sigobj.minhash.copy_and_clear()

    notify(f"loaded signature from '{args.query}'")
    notify(f"ksize={query_ksize} moltype={sigobj.minhash.moltype}")

    # track found kmers
    found_kmers = {}

    # now, iterate over the input sequences and output those that overlap
    # with hashes!
    n_seq = 0
    n = 0 # bp loaded
    m = 0 # bp in found sequences
    p = 0 # number of k-mers found
    watermark = 0
    for filename in args.seqfiles:
        for record in screed.open(filename):
            found = False
            for name, sequence in iterate_sequence_chunks(record, query_ksize):
                n += len(sequence)
                n_seq += 1
                while n >= watermark:
                    notify(f'... {filename[:20]} read {int(watermark/1e3)} kb / found kmers:{p}/{len(sigobj.minhash)}', end='\r')
                    watermark += NOTIFY_EVERY_BP

                # now do the hard work of finding the matching k-mers!
                new_mh.add_sequence(sequence, force=True)
                for kmer, hashval in get_kmers_for_hashvals(sequence,
                                                            query_hashvals,
                                                            query_ksize):
                    found_kmers[kmer] = hashval
                    found = True

                p = len(found_kmers)

            if found:
                # write out sequence?
                if seqout_fp:
                    seqout_fp.write('>{}\n{}\n'.format(name, record.sequence))
                    m += len(sequence)

    notify('done reading sequences!')

    if seqout_fp:
        notify('read {} bp, wrote {} bp in matching sequences', n, m)

    if kmerout_fp:
        for kmer, hashval in found_kmers.items():
            kmerout_w.writerow([kmer, str(hashval)])
        notify('read {} bp, found {} kmers matching hashvals', n,
               len(found_kmers))

    notify(f'jaccard similarity between query and found k-mers is {new_mh.jaccard(sigobj.minhash)}')
    notify(f'jaccard containment of found k-mers by query is {sigobj.minhash.contained_by(new_mh)}')
    notify(f'jaccard containment of query by found k-mers is {new_mh.contained_by(sigobj.minhash)}')



if __name__ == '__main__':
    sys.exit(main())
