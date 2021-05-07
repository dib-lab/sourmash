#! /usr/bin/env python3
"""
Given a list of hash values, create a sourmash signature.

CTB: Should eventually be added to sourmash signature import/export.
CTB: recommend people use scaled=1; make default?
"""
import sys
import argparse
import sourmash
from sourmash import MinHash, SourmashSignature
from sourmash import sourmash_args
from sourmash.logging import notify, error


def main():
    p = argparse.ArgumentParser()
    p.add_argument('hashfile') 					# file that contains hashes
    p.add_argument('-o', '--output', default=None,
                   help='file to output signature to')
    p.add_argument('-k', '--ksize', default=None, type=int)
    p.add_argument('--scaled', default=None, type=int)
    p.add_argument('--num', default=None, type=int)
    p.add_argument('--name', default='', help='signature name')
    p.add_argument('--filename', default='',
                   help='filename to add to signature')
    args = p.parse_args()

    # check arguments.
    if args.scaled and args.num:
        error('cannot specify both --num and --scaled! exiting.')
        return -1

    if not args.ksize:
        error('must specify --ksize')
        return -1

    if not args.output:
        error('must specify --output')
        return -1

    # first, load in all the hashes
    hashes = set()
    for line in open(args.hashfile, 'rt'):
        hashval = int(line.strip())
        hashes.add(hashval)

    if not hashes:
        error("ERROR, no hashes loaded from {}!", args.hashfile)
        return -1

    notify('loaded {} distinct hashes from {}', len(hashes), args.hashfile)

    # now, create the MinHash object that we'll use.
    scaled = 0
    num = 0
    if args.scaled:
        scaled = args.scaled
    elif args.num:
        num = args.num
    else:
        notify('setting --num automatically from the number of hashes.')
        num = len(hashes)

    # construct empty MinHash object according to args
    minhash = MinHash(n=num, ksize=args.ksize, scaled=scaled)

    # add hashes into!
    minhash.add_many(hashes)

    if len(minhash) < len(hashes):
        notify("WARNING: loaded {} hashes, but only {} made it into MinHash.",
               len(hashes), len(minhash))
        if scaled:
            notify("This is probably because of the scaled argument.")
        elif args.num:
            notify("This is probably because your --num is set to {}",
                   args.num)

    if num > len(minhash):
        notify("WARNING: --num set to {}, but only {} hashes in signature.",
               num, len(minhash))

    sigobj = sourmash.SourmashSignature(minhash, name=args.name,
                                        filename=args.filename)

    with open(args.output, 'wt') as fp:
        sourmash.save_signatures([sigobj], fp)
    notify('wrote signature to {}', args.output)


if __name__ == '__main__':
    sys.exit(main())
