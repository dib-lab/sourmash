"""
Command-line entry point for 'python -m sourmash.sig'
"""
import sys
import csv
import json
import os
from collections import defaultdict

import sourmash
import copy
from sourmash.sourmash_args import FileOutput

from sourmash.logging import set_quiet, error, notify, set_quiet, print_results, debug
from sourmash import sourmash_args
from sourmash.minhash import _get_max_hash_for_scaled

usage='''
sourmash signature <command> [<args>] - manipulate/work with signature files.

** Commands can be:

cat <signature> [<signature> ... ]        - concatenate all signatures
describe <signature> [<signature> ... ]   - show details of signature
downsample <signature> [<signature> ... ] - downsample one or more signatures
extract <signature> [<signature> ... ]    - extract one or more signatures
filter <signature> [<signature> ... ]     - filter k-mers on abundance
flatten <signature> [<signature> ... ]    - remove abundances
intersect <signature> [<signature> ...]   - intersect one or more signatures
merge <signature> [<signature> ...]       - merge one or more signatures
rename <signature> <name>                 - rename signature
split <signatures> [<signature> ...]      - split signatures into single files
subtract <signature> <other_sig> [...]    - subtract one or more signatures
import [ ... ]                            - import a mash or other signature
export <signature>                        - export a signature, e.g. to mash
overlap <signature1> <signature2>         - see detailed comparison of sigs
abundhist <signature> [<sig> ...]         - calculate abundance histogram

** Use '-h' to get subcommand-specific help, e.g.

sourmash signature merge -h
'''


def _check_abundance_compatibility(sig1, sig2):
    if sig1.minhash.track_abundance != sig2.minhash.track_abundance:
        raise ValueError("incompatible signatures: track_abundance is {} in first sig, {} in second".format(sig1.minhash.track_abundance, sig2.minhash.track_abundance))


def _set_num_scaled(mh, num, scaled):
    "set num and scaled values on a MinHash object"
    mh_params = list(mh.__getstate__())
    # Number of hashes is 0th parameter
    mh_params[0] = num
    # Scale is 8th parameter
    mh_params[8] = _get_max_hash_for_scaled(scaled)
    mh.__setstate__(mh_params)
    assert mh.num == num
    assert mh.scaled == scaled


##### actual command line functions


def cat(args):
    """
    concatenate all signatures into one file.
    """
    set_quiet(args.quiet)

    encountered_md5sums = defaultdict(int)   # used by --unique
    progress = sourmash_args.SignatureLoadingProgress()

    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()

    for sigfile in args.signatures:
        try:
            loader = sourmash_args.load_file_as_signatures(sigfile,
                                                           progress=progress)
            n_loaded = 0
            for sig in loader:
                n_loaded += 1

                md5 = sig.md5sum()
                encountered_md5sums[md5] += 1
                if args.unique and encountered_md5sums[md5] > 1:
                    continue

                save_sigs.add(sig)
        except Exception as exc:
            error(str(exc))
            error('(continuing)')

        notify('loaded {} signatures from {}...', n_loaded, sigfile, end='\r')

    notify('loaded {} signatures total.', len(save_sigs))

    save_sigs.close()

    notify('output {} signatures', len(save_sigs))

    multiple_md5 = [ 1 for cnt in encountered_md5sums.values() if cnt > 1 ]
    if multiple_md5:
        notify('encountered {} MinHashes multiple times', sum(multiple_md5))
        if args.unique:
            notify('...and removed the duplicates, because --unique was specified.')


def split(args):
    """
    split all signatures into individual files
    """
    set_quiet(args.quiet)

    output_names = set()
    output_scaled_template = '{md5sum}.k={ksize}.scaled={scaled}.{moltype}.dup={dup}.{basename}.sig'
    output_num_template = '{md5sum}.k={ksize}.num={num}.{moltype}.dup={dup}.{basename}.sig'

    if args.outdir:
        if not os.path.exists(args.outdir):
            notify('Creating --outdir {}', args.outdir)
            os.mkdir(args.outdir)

    progress = sourmash_args.SignatureLoadingProgress()

    total = 0
    for sigfile in args.signatures:
        # load signatures from input file:
        this_siglist = sourmash_args.load_file_as_signatures(sigfile,
                                                             progress=progress)

        # save each file individually --
        n_signatures = 0
        for sig in this_siglist:
            n_signatures += 1
            md5sum = sig.md5sum()[:8]
            minhash = sig.minhash
            basename = os.path.basename(sig.filename)
            if not basename or basename == '-':
                basename = 'none'

            params = dict(basename=basename,
                          md5sum=md5sum,
                          scaled=minhash.scaled,
                          ksize=minhash.ksize,
                          num=minhash.num,
                          moltype=minhash.moltype)

            if minhash.scaled:
                output_template = output_scaled_template
            else: # num
                assert minhash.num
                output_template = output_num_template

            # figure out if this is duplicate, build unique filename
            n = 0
            params['dup'] = n
            output_name = output_template.format(**params)
            while output_name in output_names:
                params['dup'] = n
                output_name = output_template.format(**params)
                n += 1

            output_names.add(output_name)

            if args.outdir:
                output_name = os.path.join(args.outdir, output_name)

            if os.path.exists(output_name):
                notify("** overwriting existing file {}".format(output_name))

            # save!
            with open(output_name, 'wt') as outfp:
                sourmash.save_signatures([sig], outfp)
                notify('writing sig to {}', output_name)

        notify('loaded {} signatures from {}...', n_signatures, sigfile,
               end='\r')
        total += n_signatures

    notify('loaded and split {} signatures total.', total)


def describe(args):
    """
    provide basic info on signatures
    """
    set_quiet(args.quiet)

    # write CSV?
    w = None
    csv_fp = None
    if args.csv:
        csv_fp = open(args.csv, 'w', newline='')
        w = csv.DictWriter(csv_fp,
                           ['signature_file', 'md5', 'ksize', 'moltype', 'num',
                            'scaled', 'n_hashes', 'seed', 'with_abundance',
                            'name', 'filename', 'license'],
                           extrasaction='ignore')
        w.writeheader()

    # load signatures and display info.
    progress = sourmash_args.SignatureLoadingProgress()

    n_loaded = 0
    for signature_file in args.signatures:
        try:
            loader = sourmash_args.load_file_as_signatures(signature_file,
                                                           progress=progress)
            for sig in loader:
                n_loaded += 1

                # extract info, write as appropriate.
                mh = sig.minhash
                ksize = mh.ksize
                moltype = mh.moltype
                scaled = mh.scaled
                num = mh.num
                seed = mh.seed
                n_hashes = len(mh)
                with_abundance = 0
                if mh.track_abundance:
                    with_abundance = 1
                md5 = sig.md5sum()
                name = sig.name or "** no name **"
                filename = sig.filename or "** no name **"
                license = sig.license

                if w:
                    w.writerow(locals())

                print_results('''\
---
signature filename: {signature_file}
signature: {name}
source file: {filename}
md5: {md5}
k={ksize} molecule={moltype} num={num} scaled={scaled} seed={seed} track_abundance={with_abundance}
size: {n_hashes}
signature license: {license}
''', **locals())

        except Exception as exc:
            error('\nError while reading signatures from {}:'.format(signature_file))
            error(str(exc))
            error('(continuing)')
            raise

    notify('loaded {} signatures total.', n_loaded)

    if csv_fp:
        csv_fp.close()


def overlap(args):
    """
    provide detailed comparison of two signatures
    """
    set_quiet(args.quiet)

    moltype = sourmash_args.calculate_moltype(args)

    sig1 = sourmash.load_one_signature(args.signature1, ksize=args.ksize,
                                       select_moltype=moltype)
    sig2 = sourmash.load_one_signature(args.signature2, ksize=args.ksize,
                                       select_moltype=moltype)

    notify('loaded one signature each from {} and {}', args.signature1,
           args.signature2)

    try:
        similarity = sig1.similarity(sig2)
    except ValueError:
        raise

    cont1 = sig1.contained_by(sig2)
    cont2 = sig2.contained_by(sig1)

    sig1_file = args.signature1
    sig2_file = args.signature2

    name1 = sig1.name
    name2 = sig2.name

    md5_1 = sig1.md5sum()
    md5_2 = sig2.md5sum()

    ksize = sig1.minhash.ksize
    moltype = sig1.minhash.moltype

    num = sig1.minhash.num
    size1 = len(sig1.minhash)
    size2 = len(sig2.minhash)

    scaled = sig1.minhash.scaled

    hashes_1 = set(sig1.minhash.hashes)
    hashes_2 = set(sig2.minhash.hashes)

    num_common = len(hashes_1.intersection(hashes_2))
    disjoint_1 = len(hashes_1 - hashes_2)
    disjoint_2 = len(hashes_2 - hashes_1)
    num_union = len(hashes_1.union(hashes_2))

    print('''\
first signature:
  signature filename: {sig1_file}
  signature: {name1}
  md5: {md5_1}
  k={ksize} molecule={moltype} num={num} scaled={scaled}

second signature:
  signature filename: {sig2_file}
  signature: {name2}
  md5: {md5_2}
  k={ksize} molecule={moltype} num={num} scaled={scaled}

similarity:                  {similarity:.5f}
first contained in second:   {cont1:.5f}
second contained in first:   {cont2:.5f}

number of hashes in first:   {size1}
number of hashes in second:  {size2}

number of hashes in common:  {num_common}
only in first:               {disjoint_1}
only in second:              {disjoint_2}
total (union):               {num_union}
'''.format(**locals()))


def merge(args):
    """
    merge one or more signatures.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    first_sig = None
    mh = None
    total_loaded = 0

    # iterate over all the sigs from all the files.
    progress = sourmash_args.SignatureLoadingProgress()

    for sigfile in args.signatures:
        notify('loading signatures from {}...', sigfile, end='\r')
        this_n = 0
        for sigobj in sourmash_args.load_file_as_signatures(sigfile,
                                                        ksize=args.ksize,
                                                        select_moltype=moltype,
                                                        progress=progress):

            # first signature? initialize a bunch of stuff
            if first_sig is None:
                first_sig = sigobj
                mh = first_sig.minhash.copy_and_clear()

                # forcibly remove abundance?
                if args.flatten:
                    mh.track_abundance = False

            try:
                sigobj_mh = sigobj.minhash
                if not args.flatten:
                    _check_abundance_compatibility(first_sig, sigobj)
                else:
                    sigobj_mh.track_abundance = False

                mh.merge(sigobj_mh)
            except:
                error("ERROR when merging signature '{}' ({}) from file {}",
                      sigobj, sigobj.md5sum()[:8], sigfile)
                raise

            this_n += 1
            total_loaded += 1
        if this_n:
            notify('loaded and merged {} signatures from {}...', this_n, sigfile, end='\r')

    if not total_loaded:
        error("no signatures to merge!?")
        sys.exit(-1)

    merged_sigobj = sourmash.SourmashSignature(mh, name=args.name)

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures([merged_sigobj], fp=fp)

    notify('loaded and merged {} signatures', total_loaded)


def intersect(args):
    """
    intersect one or more signatures by taking the intersection of hashes.

    This function always removes abundances.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    first_sig = None
    mins = None
    total_loaded = 0

    progress = sourmash_args.SignatureLoadingProgress()

    for sigfile in args.signatures:
        for sigobj in sourmash_args.load_file_as_signatures(sigfile,
                                               ksize=args.ksize,
                                               select_moltype=moltype,
                                               progress=progress):
            if first_sig is None:
                first_sig = sigobj
                mins = set(sigobj.minhash.hashes)
            else:
                # check signature compatibility --
                if not sigobj.minhash.is_compatible(first_sig.minhash):
                    error("incompatible minhashes; specify -k and/or molecule type.")
                    sys.exit(-1)

            mins.intersection_update(sigobj.minhash.hashes)
            total_loaded += 1
        notify('loaded and intersected signatures from {}...', sigfile, end='\r')

    if total_loaded == 0:
        error("no signatures to merge!?")
        sys.exit(-1)

    # forcibly turn off track_abundance, unless --abundances-from set.
    if not args.abundances_from:
        intersect_mh = first_sig.minhash.copy_and_clear()
        intersect_mh.track_abundance = False
        intersect_mh.add_many(mins)
        intersect_sigobj = sourmash.SourmashSignature(intersect_mh)
    else:
        notify('loading signature from {}, keeping abundances',
               args.abundances_from)
        abund_sig = sourmash.load_one_signature(args.abundances_from,
                                                ksize=args.ksize,
                                                select_moltype=moltype)
        if not abund_sig.minhash.track_abundance:
            error("--track-abundance not set on loaded signature?! exiting.")
            sys.exit(-1)
        intersect_mh = abund_sig.minhash.copy_and_clear()
        abund_mins = abund_sig.minhash.hashes

        # do one last intersection
        mins.intersection_update(abund_mins)
        abund_mins = { k: abund_mins[k] for k in mins }

        intersect_mh.set_abundances(abund_mins)
        intersect_sigobj = sourmash.SourmashSignature(intersect_mh)

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures([intersect_sigobj], fp=fp)

    notify('loaded and intersected {} signatures', total_loaded)


def subtract(args):
    """
    subtract one or more signatures from another
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    from_sigfile = args.signature_from
    from_sigobj = sourmash.load_one_signature(from_sigfile, ksize=args.ksize, select_moltype=moltype)

    from_mh = from_sigobj.minhash
    if from_mh.track_abundance and not args.flatten:
        error('Cannot use subtract on signatures with abundance tracking, sorry!')
        sys.exit(1)

    subtract_mins = set(from_mh.hashes)

    notify('loaded signature from {}...', from_sigfile, end='\r')

    progress = sourmash_args.SignatureLoadingProgress()

    total_loaded = 0
    for sigfile in args.subtraction_sigs:
        for sigobj in sourmash_args.load_file_as_signatures(sigfile,
                                                        ksize=args.ksize,
                                                        select_moltype=moltype,
                                                        progress=progress):
            if not sigobj.minhash.is_compatible(from_mh):
                error("incompatible minhashes; specify -k and/or molecule type.")
                sys.exit(-1)

            if sigobj.minhash.track_abundance and not args.flatten:
                error('Cannot use subtract on signatures with abundance tracking, sorry!')
                sys.exit(1)

            subtract_mins -= set(sigobj.minhash.hashes)

            notify('loaded and subtracted signatures from {}...', sigfile, end='\r')
            total_loaded += 1

    if not total_loaded:
        error("no signatures to subtract!?")
        sys.exit(-1)


    subtract_mh = from_sigobj.minhash.copy_and_clear()
    subtract_mh.add_many(subtract_mins)

    subtract_sigobj = sourmash.SourmashSignature(subtract_mh)

    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures([subtract_sigobj], fp=fp)

    notify('loaded and subtracted {} signatures', total_loaded)


def rename(args):
    """
    rename one or more signatures.
    """
    set_quiet(args.quiet, args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    progress = sourmash_args.SignatureLoadingProgress()

    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()

    for filename in args.sigfiles:
        debug('loading {}', filename)
        siglist = sourmash_args.load_file_as_signatures(filename,
                                                        ksize=args.ksize,
                                                        select_moltype=moltype,
                                                        progress=progress)

        for sigobj in siglist:
            sigobj._name = args.name
            save_sigs.add(sigobj)

    save_sigs.close()

    notify("set name to '{}' on {} signatures", args.name, len(save_sigs))


def extract(args):
    """
    extract signatures.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    progress = sourmash_args.SignatureLoadingProgress()

    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()

    total_loaded = 0
    for filename in args.signatures:
        siglist = sourmash_args.load_file_as_signatures(filename,
                                                        ksize=args.ksize,
                                                        select_moltype=moltype,
                                                        progress=progress)
        siglist = list(siglist)

        total_loaded += len(siglist)

        # select!
        if args.md5 is not None:
            siglist = [ ss for ss in siglist if args.md5 in ss.md5sum() ]
        if args.name is not None:
            siglist = [ ss for ss in siglist if args.name in str(ss) ]

        for ss in siglist:
            save_sigs.add(ss)

    notify("loaded {} total that matched ksize & molecule type",
           total_loaded)
    if not save_sigs:
        error("no matching signatures!")
        sys.exit(-1)

    save_sigs.close()

    notify("extracted {} signatures from {} file(s)", len(save_sigs),
           len(args.signatures))


def filter(args):
    """
    filter hashes by abundance in all of the signatures
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    progress = sourmash_args.SignatureLoadingProgress()

    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()

    total_loaded = 0
    for filename in args.signatures:
        siglist = sourmash_args.load_file_as_signatures(filename,
                                                        ksize=args.ksize,
                                                        select_moltype=moltype,
                                                        progress=progress)
        siglist = list(siglist)

        total_loaded += len(siglist)

        # select!
        if args.md5 is not None:
            siglist = [ ss for ss in siglist if args.md5 in ss.md5sum() ]
        if args.name is not None:
            siglist = [ ss for ss in siglist if args.name in str(ss) ]

        for ss in siglist:
            mh = ss.minhash
            if not mh.track_abundance:
                notify('ignoring signature {} - track_abundance not set.',
                       ss)
                continue

            abunds = mh.hashes
            abunds2 = {}
            for k, v in abunds.items():
                if v >= args.min_abundance:
                    if args.max_abundance is None or \
                       v <= args.max_abundance:
                       abunds2[k] = v

            filtered_mh = mh.copy_and_clear()
            filtered_mh.set_abundances(abunds2)

            ss.minhash = filtered_mh

            save_sigs.add(ss)

    save_sigs.close()

    notify("loaded {} total that matched ksize & molecule type",
           total_loaded)
    notify("extracted {} signatures from {} file(s)", len(save_sigs),
           len(args.signatures))


def flatten(args):
    """
    flatten one or more signatures, removing abundances.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    progress = sourmash_args.SignatureLoadingProgress()

    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()

    total_loaded = 0
    for filename in args.signatures:
        siglist = sourmash_args.load_file_as_signatures(filename,
                                                        ksize=args.ksize,
                                                        select_moltype=moltype,
                                                        progress=progress)
        siglist = list(siglist)

        total_loaded += len(siglist)

        # select!
        if args.md5 is not None:
            siglist = [ ss for ss in siglist if args.md5 in ss.md5sum() ]
        if args.name is not None:
            siglist = [ ss for ss in siglist if args.name in ss.name ]

        for ss in siglist:
            ss.minhash = ss.minhash.flatten()
            save_sigs.add(ss)

    save_sigs.close()

    notify("loaded {} total that matched ksize & molecule type",
           total_loaded)
    notify("extracted {} signatures from {} file(s)", len(save_sigs),
           len(args.signatures))


def downsample(args):
    """
    downsample a scaled signature.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    if not args.num and not args.scaled:
        error('must specify either --num or --scaled value')
        sys.exit(-1)

    if args.num and args.scaled:
        error('cannot specify both --num and --scaled')
        sys.exit(-1)

    save_sigs = sourmash_args.SaveSignaturesToLocation(args.output)
    save_sigs.open()

    progress = sourmash_args.SignatureLoadingProgress()

    total_loaded = 0
    for sigfile in args.signatures:
        siglist = sourmash_args.load_file_as_signatures(sigfile,
                                                        ksize=args.ksize,
                                                        select_moltype=moltype,
                                                        progress=progress)

        for sigobj in siglist:
            mh = sigobj.minhash

            notify('loading and downsampling signature from {}...', sigfile, end='\r')
            total_loaded += 1
            if args.scaled:
                if mh.scaled:
                    mh_new = mh.downsample(scaled=args.scaled)
                else:                         # try to turn a num into a scaled
                    # first check: can we?
                    max_hash = _get_max_hash_for_scaled(args.scaled)
                    mins = mh.hashes
                    if max(mins) < max_hash:
                        raise ValueError("this num MinHash does not have enough hashes to convert it into a scaled MinHash.")

                    mh_new = copy.copy(mh)
                    _set_num_scaled(mh_new, 0, args.scaled)
            elif args.num:
                if mh.num:
                    mh_new = mh.downsample(num=args.num)
                else:                         # try to turn a scaled into a num
                    # first check: can we?
                    if len(mh) < args.num:
                        raise ValueError("this scaled MinHash has only {} hashes")

                    mh_new = copy.copy(mh)
                    _set_num_scaled(mh_new, args.num, 0)

            sigobj.minhash = mh_new

            save_sigs.add(sigobj)

    save_sigs.close()

    notify("loaded and downsampled {} signatures", total_loaded)


def sig_import(args):
    """
    import a signature into sourmash format.
    """
    set_quiet(args.quiet)

    siglist = []
    if args.csv:
        for filename in args.filenames:
            with open(filename, newline='') as csv_fp:
                reader = csv.reader(csv_fp)
                siglist = []
                for row in reader:
                    hashfn = row[0]
                    hashseed = int(row[1])

                    # only support a limited import type, for now ;)
                    assert hashfn == 'murmur64'
                    assert hashseed == 42

                    _, _, ksize, name, hashes = row
                    ksize = int(ksize)

                    hashes = hashes.strip()
                    hashes = list(map(int, hashes.split(' ' )))

                    e = sourmash.MinHash(len(hashes), ksize)
                    e.add_many(hashes)
                    s = sourmash.SourmashSignature(e, filename=name)
                    siglist.append(s)
                    notify('loaded signature: {} {}', name, s.md5sum()[:8])
    else:
        for filename in args.filenames:
            with open(filename) as fp:
                x = json.loads(fp.read())

            ksize = x['kmer']
            num = x['sketchSize']

            assert x['hashType'] == "MurmurHash3_x64_128"
            assert x['hashBits'] == 64
            assert x['hashSeed'] == 42

            xx = x['sketches'][0]
            hashes = xx['hashes']

            mh = sourmash.MinHash(ksize=ksize, n=num, is_protein=False)
            mh.add_many(hashes)

            s = sourmash.SourmashSignature(mh, filename=filename)
            siglist.append(s)

    notify('saving {} signatures to JSON', len(siglist))
    with FileOutput(args.output, 'wt') as fp:
        sourmash.save_signatures(siglist, fp)


def export(args):
    """
    export a signature to mash format
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    query = sourmash_args.load_query_signature(args.filename,
                                               ksize=args.ksize,
                                               select_moltype=moltype,
                                               select_md5=args.md5)
    mh = query.minhash

    x = {}
    x['kmer'] = mh.ksize
    x['sketchSize'] = len(mh)

    x['hashType'] = "MurmurHash3_x64_128"
    x['hashBits'] = 64
    x['hashSeed'] = mh.seed

    ll = list(mh.hashes)
    x['sketches'] = [{ 'hashes': ll }]

    with FileOutput(args.output, 'wt') as fp:
        print(json.dumps(x), file=fp)
    notify("exported signature {} ({})", query, query.md5sum()[:8])




def abundhist(args):
    """
    output abundance histogram and/or raw abundances.
    """
    import numpy, collections
    import termplotlib as tpl

    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    outlist = []
    total_loaded = 0
    for filename in args.signatures:
        siglist = sourmash.load_file_as_signatures(filename, ksize=args.ksize,
                                                   select_moltype=moltype)
        siglist = list(siglist)

        total_loaded += len(siglist)

        # select!
        if args.md5 is not None:
            siglist = [ ss for ss in siglist if args.md5 in ss.md5sum() ]
        if args.name is not None:
            siglist = [ ss for ss in siglist if args.name in ss.name() ]

    notify("loaded {} total that matched ksize & molecule type",
           total_loaded)
    if len(siglist) != total_loaded:
        notify("selected {} via name / md5 selectors".format(len(siglist)))
    notify('')

    counts_d = collections.defaultdict(int)
    for ss in siglist:
        for hashval, abund in ss.minhash.hashes.items():
            counts_d[hashval] += abund

    all_counts = list(counts_d.values())

    min_range = 1
    if args.min is not None:
        min_range = args.min
    max_range = max(all_counts)
    if args.max is not None:
        max_range = args.max

    n_bins = args.bins
    if max_range - min_range + 1 < n_bins:
        n_bins = max_range - min_range + 1

    # make hist
    counts, bin_edges = numpy.histogram(all_counts,
                                        range=(min_range, max_range),
                                        bins=n_bins)
    bin_edges = bin_edges.astype(int)

    # plot
    fig = tpl.figure()
    f = fig.barh(counts, [ str(x) for x in bin_edges[1:] ], force_ascii=True)
    fig.show()

    # output histogram in csv?
    if args.output:
        with FileOutput(args.output, 'wt') as fp:
            w = csv.writer(fp)
            w.writerow(['count', 'n_count'])
            for nc, c in zip(counts, bin_edges[1:]):
                w.writerow([c, nc])

    # output raw counts tagged with hashval?
    if args.abundances:
        with FileOutput(args.abundances, 'wt') as fp:
            w = csv.writer(fp)
            w.writerow(['hashval', 'count'])
            for hashval, count in counts_d.items():
                w.writerow([hashval, count])


def main(arglist=None):
    args = sourmash.cli.get_parser().parse_args(arglist)
    submod = getattr(sourmash.cli.sig, args.subcmd)
    mainmethod = getattr(submod, 'main')
    return mainmethod(args)


if __name__ == '__main__':
    main(sys.argv)
