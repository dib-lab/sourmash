"""
Functions implementing the main command-line subcommands.
"""
from __future__ import print_function, division, absolute_import

import csv
import os
import os.path
import sys

import screed
from .compare import compare_all_pairs, compare_serial_containment
from . import MinHash, load_sbt_index, create_sbt_index
from . import signature as sig
from . import sourmash_args
from .logging import notify, error, print_results, set_quiet
from .sbtmh import SearchMinHashesFindBest, SigLeaf

from .sourmash_args import DEFAULT_LOAD_K, FileOutput

DEFAULT_N = 500
WATERMARK_SIZE = 10000

from .command_compute import compute


def compare(args):
    "Compare multiple signature files and create a distance matrix."
    import numpy

    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    # check directories for all signatures
    if args.traverse_directory:
        inp_files = list(sourmash_args.traverse_find_sigs(args.signatures))
    else:
        inp_files = list(args.signatures)

    # load in the various signatures
    siglist = []
    ksizes = set()
    moltypes = set()
    for filename in inp_files:
        if not os.path.exists(filename) and not \
               (args.force or args.traverse_directory):
            error("file '{}' does not exist! exiting.", filename)
            sys.exit(-1)

        notify("loading '{}'", filename, end='\r')
        loaded = sig.load_signatures(filename,
                                     ksize=args.ksize,
                                     select_moltype=moltype)
        loaded = list(loaded)
        if not loaded:
            notify('\nwarning: no signatures loaded at given ksize/molecule type from {}', filename)
        siglist.extend(loaded)

        # track ksizes/moltypes
        for s in loaded:
            ksizes.add(s.minhash.ksize)
            moltypes.add(sourmash_args.get_moltype(s))

        # error out while loading if we have more than one ksize/moltype
        if len(ksizes) > 1 or len(moltypes) > 1:
            break

    if not siglist:
        error('no signatures found! exiting.')
        sys.exit(-1)

    # check ksizes and type
    if len(ksizes) > 1:
        error('multiple k-mer sizes loaded; please specify one with -k.')
        ksizes = sorted(ksizes)
        error('(saw k-mer sizes {})'.format(', '.join(map(str, ksizes))))
        sys.exit(-1)

    if len(moltypes) > 1:
        error('multiple molecule types loaded; please specify --dna, --protein')
        sys.exit(-1)

    notify(' '*79, end='\r')
    notify('loaded {} signatures total.'.format(len(siglist)))

    # check to make sure they're potentially compatible - either using
    # max_hash/scaled, or not.
    scaled_sigs = [s.minhash.max_hash for s in siglist]
    is_scaled = all(scaled_sigs)
    is_scaled_2 = any(scaled_sigs)

    # complain if it's not all one or the other
    if is_scaled != is_scaled_2:
        error('cannot mix scaled signatures with bounded signatures')
        sys.exit(-1)

    # complain if --containment and not is_scaled
    if args.containment and not is_scaled:
        error('must use scaled signatures with --containment option')
        sys.exit(-1)

    # notify about implicit --ignore-abundance:
    if args.containment:
        track_abundances = any(( s.minhash.track_abundance for s in siglist ))
        if track_abundances:
            notify('NOTE: --containment means signature abundances are flattened.')

    # if using --scaled, downsample appropriately
    printed_scaled_msg = False
    if is_scaled:
        max_scaled = max(s.minhash.scaled for s in siglist)
        for s in siglist:
            if s.minhash.scaled != max_scaled:
                if not printed_scaled_msg:
                    notify('downsampling to scaled value of {}'.format(max_scaled))
                    printed_scaled_msg = True
                s.minhash = s.minhash.downsample_scaled(max_scaled)

    if len(siglist) == 0:
        error('no signatures!')
        sys.exit(-1)

    notify('')

    # build the distance matrix
    numpy.set_printoptions(precision=3, suppress=True)

    # do all-by-all calculation

    labeltext = [item.name() for item in siglist]
    if args.containment:
        similarity = compare_serial_containment(siglist)
    else:
        similarity = compare_all_pairs(siglist, args.ignore_abundance,
                                       n_jobs=args.processes)

    if len(siglist) < 30:
        for i, E in enumerate(siglist):
            # for small matrices, pretty-print some output
            name_num = '{}-{}'.format(i, E.name())
            if len(name_num) > 20:
                name_num = name_num[:17] + '...'
            print_results('{:20s}\t{}'.format(name_num, similarity[i, :, ],))

    print_results('min similarity in matrix: {:.3f}', numpy.min(similarity))
    # shall we output a matrix?
    if args.output:
        labeloutname = args.output + '.labels.txt'
        notify('saving labels to: {}', labeloutname)
        with open(labeloutname, 'w') as fp:
            fp.write("\n".join(labeltext))

        notify('saving comparison matrix to: {}', args.output)
        with open(args.output, 'wb') as fp:
            numpy.save(fp, similarity)

    # output CSV?
    if args.csv:
        with FileOutput(args.csv, 'wt') as csv_fp:
            w = csv.writer(csv_fp)
            w.writerow(labeltext)

            for i in range(len(labeltext)):
                y = []
                for j in range(len(labeltext)):
                    y.append('{}'.format(similarity[i][j]))
                w.writerow(y)


def plot(args):
    "Produce a clustering matrix and plot."
    import matplotlib as mpl
    mpl.use('Agg')
    import numpy
    import pylab
    import scipy.cluster.hierarchy as sch
    from . import fig as sourmash_fig

    # load files
    D_filename = args.distances
    labelfilename = D_filename + '.labels.txt'

    notify('loading comparison matrix from {}...', D_filename)
    D = numpy.load(open(D_filename, 'rb'))
    notify('...got {} x {} matrix.', *D.shape)

    if args.labeltext:
        labelfilename = args.labeltext
    notify('loading labels from {}', labelfilename)
    labeltext = [ x.strip() for x in open(labelfilename) ]
    if len(labeltext) != D.shape[0]:
        error('{} labels != matrix size, exiting', len(labeltext))
        sys.exit(-1)

    # build filenames, decide on PDF/PNG output
    dendrogram_out = os.path.basename(D_filename) + '.dendro'
    if args.pdf:
        dendrogram_out += '.pdf'
    else:
        dendrogram_out += '.png'

    matrix_out = os.path.basename(D_filename) + '.matrix'
    if args.pdf:
        matrix_out += '.pdf'
    else:
        matrix_out += '.png'

    hist_out = os.path.basename(D_filename) + '.hist'
    if args.pdf:
        hist_out += '.pdf'
    else:
        hist_out += '.png'

    # output to a different directory?
    if args.output_dir:
        if not os.path.isdir(args.output_dir):
            os.mkdir(args.output_dir)
        dendrogram_out = os.path.join(args.output_dir, dendrogram_out)
        matrix_out = os.path.join(args.output_dir, matrix_out)
        hist_out = os.path.join(args.output_dir, hist_out)

    # make the histogram
    notify('saving histogram of matrix values => {}', hist_out)
    fig = pylab.figure(figsize=(8,5))
    pylab.hist(numpy.array(D.flat), bins=100)
    fig.savefig(hist_out)

    ### make the dendrogram:
    fig = pylab.figure(figsize=(8,5))
    ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
    ax1.set_xticks([])
    ax1.set_yticks([])

    # subsample?
    if args.subsample:
        numpy.random.seed(args.subsample_seed)

        sample_idx = list(range(len(labeltext)))
        numpy.random.shuffle(sample_idx)
        sample_idx = sample_idx[:args.subsample]

        np_idx = numpy.array(sample_idx)
        D = D[numpy.ix_(np_idx, np_idx)]
        labeltext = [ labeltext[idx] for idx in sample_idx ]

    ### do clustering
    Y = sch.linkage(D, method='single')
    sch.dendrogram(Y, orientation='right', labels=labeltext)
    fig.savefig(dendrogram_out)
    notify('wrote dendrogram to: {}', dendrogram_out)

    ### make the dendrogram+matrix:
    fig = sourmash_fig.plot_composite_matrix(D, labeltext,
                                             show_labels=args.labels,
                                             show_indices=args.indices,
                                             vmin=args.vmin,
                                             vmax=args.vmax,
                                             force=args.force)
    fig.savefig(matrix_out)
    notify('wrote numpy distance matrix to: {}', matrix_out)

    if len(labeltext) < 30:
        # for small matrices, print out sample numbering for FYI.
        for i, name in enumerate(labeltext):
            print_results('{}\t{}', i, name)


def import_csv(args):
    "Import a CSV file full of signatures/hashes."

    with open(args.mash_csvfile, 'r') as fp:
        reader = csv.reader(fp)
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

            e = MinHash(len(hashes), ksize)
            e.add_many(hashes)
            s = sig.SourmashSignature(e, filename=name)
            siglist.append(s)
            notify('loaded signature: {} {}', name, s.md5sum()[:8])

        notify('saving {} signatures to JSON', len(siglist))
        with FileOutput(args.output, 'wt') as outfp:
            sig.save_signatures(siglist, outfp)


def dump(args):
    "Dump hashes for each input signature into a {name}.dump.txt file."
    for filename in args.filenames:
        notify('loading {}', filename)
        siglist = sig.load_signatures(filename, ksize=args.ksize)
        siglist = list(siglist)
        assert len(siglist) == 1

        s = siglist[0]

        fp = open(filename + '.dump.txt', 'w')
        fp.write(" ".join((map(str, s.minhash.get_hashes()))))
        fp.close()


def sbt_combine(args):
    inp_files = list(args.sbts)
    notify('combining {} SBTs', len(inp_files))

    tree = load_sbt_index(inp_files.pop(0))

    for f in inp_files:
        new_tree = load_sbt_index(f)
        # TODO: check if parameters are the same for both trees!
        tree.combine(new_tree)

    notify('saving SBT under "{}".', args.sbt_name)
    tree.save(args.sbt_name)


def index(args):
    """
    Build a Sequence Bloom Tree index of the given signatures.
    """
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    if args.append:
        tree = load_sbt_index(args.sbt_name)
    else:
        tree = create_sbt_index(args.bf_size, n_children=args.n_children)

    if args.traverse_directory:
        inp_files = list(sourmash_args.traverse_find_sigs(args.signatures,
                                                          args.force))
    else:
        inp_files = list(args.signatures)

    if args.sparseness < 0 or args.sparseness > 1.0:
        error('sparseness must be in range [0.0, 1.0].')

    if args.scaled:
        args.scaled = int(args.scaled)
        notify('downsampling signatures to scaled={}', args.scaled)

    notify('loading {} files into SBT', len(inp_files))

    n = 0
    ksizes = set()
    moltypes = set()
    nums = set()
    scaleds = set()
    md5sums = set()
    for f in inp_files:
        if n % 100 == 0:
            notify('\r...reading from {} ({} signatures so far)', f, n, end='')
        siglist = sig.load_signatures(f, ksize=args.ksize,
                                      select_moltype=moltype)

        # load all matching signatures in this file
        ss = None
        for ss in siglist:
            ksizes.add(ss.minhash.ksize)
            moltypes.add(sourmash_args.get_moltype(ss))
            nums.add(ss.minhash.num)

            if args.scaled:
                ss.minhash = ss.minhash.downsample_scaled(args.scaled)
            scaleds.add(ss.minhash.scaled)

            md5sum = ss.md5sum()
            if md5sum in md5sums:
                notify('leaf already exists, failing.')
                sys.exit(-1)

            tree.insert(ss)
            md5sums.add(md5sum)
            n += 1

        if not ss:
            continue

        # check to make sure we aren't loading incompatible signatures
        if len(ksizes) > 1 or len(moltypes) > 1:
            error('multiple k-mer sizes or molecule types present; fail.')
            error('specify --dna/--protein and --ksize as necessary')
            error('ksizes: {}; moltypes: {}',
                  ", ".join(map(str, ksizes)), ", ".join(moltypes))
            sys.exit(-1)

        if nums == { 0 } and len(scaleds) == 1:
            pass # good
        elif scaleds == { 0 } and len(nums) == 1:
            pass # also good
        else:
            error('trying to build an SBT with incompatible signatures.')
            error('nums = {}; scaleds = {}', repr(nums), repr(scaleds))
            sys.exit(-1)

    notify('')

    # did we load any!?
    if n == 0:
        error('no signatures found to load into tree!? failing.')
        sys.exit(-1)

    notify('loaded {} sigs; saving SBT under "{}"', n, args.sbt_name)
    tree.save(args.sbt_name, sparseness=args.sparseness)


def search(args):
    from .search import search_databases

    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    # set up the query.
    query = sourmash_args.load_query_signature(args.query,
                                               ksize=args.ksize,
                                               select_moltype=moltype)
    notify('loaded query: {}... (k={}, {})', query.name()[:30],
                                             query.minhash.ksize,
                                             sourmash_args.get_moltype(query))

    # downsample if requested
    if args.scaled:
        if query.minhash.max_hash == 0:
            error('cannot downsample a signature not created with --scaled')
            sys.exit(-1)

        if args.scaled != query.minhash.scaled:
            notify('downsampling query from scaled={} to {}',
                   query.minhash.scaled, int(args.scaled))
        query.minhash = query.minhash.downsample_scaled(args.scaled)

    # set up the search databases
    databases = sourmash_args.load_dbs_and_sigs(args.databases, query,
                                                not args.containment,
                                                args.traverse_directory)

    # forcibly ignore abundances if query has no abundances
    if not query.minhash.track_abundance:
        args.ignore_abundance = True

    if not len(databases):
        error('Nothing found to search!')
        sys.exit(-1)

    # do the actual search
    results = search_databases(query, databases,
                               args.threshold, args.containment,
                               args.best_only, args.ignore_abundance,
                               unload_data=True)

    n_matches = len(results)
    if args.best_only:
        args.num_results = 1

    if not args.num_results or n_matches <= args.num_results:
        print_results('{} matches:'.format(len(results)))
    else:
        print_results('{} matches; showing first {}:',
               len(results), args.num_results)
        n_matches = args.num_results

    # output!
    print_results("similarity   match")
    print_results("----------   -----")
    for sr in results[:n_matches]:
        pct = '{:.1f}%'.format(sr.similarity*100)
        name = sr.match._display_name(60)
        print_results('{:>6}       {}', pct, name)

    if args.best_only:
        notify("** reporting only one match because --best-only was set")

    if args.output:
        fieldnames = ['similarity', 'name', 'filename', 'md5']

        with FileOutput(args.output, 'wt') as fp:
            w = csv.DictWriter(fp, fieldnames=fieldnames)

            w.writeheader()
            for sr in results:
                d = dict(sr._asdict())
                del d['match']
                w.writerow(d)

    # save matching signatures upon request
    if args.save_matches:
        notify('saving all matched signatures to "{}"', args.save_matches)
        with FileOutput(args.save_matches, 'wt') as fp:
            sig.save_signatures([ sr.match for sr in results ], fp)


def categorize(args):
    "Use an SBT to find the best match to many signatures."
    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    already_names = set()
    if args.load_csv:
        with open(args.load_csv, 'rt') as fp:
            r = csv.reader(fp)
            for row in r:
                already_names.add(row[0])

    tree = load_sbt_index(args.sbt_name)

    if args.traverse_directory:
        inp_files = set(sourmash_args.traverse_find_sigs(args.queries))
    else:
        inp_files = set(args.queries) - already_names

    inp_files = set(inp_files) - already_names

    notify('found {} files to query', len(inp_files))

    loader = sourmash_args.LoadSingleSignatures(inp_files,
                                                args.ksize, moltype)

    csv_w = None
    csv_fp = None
    if args.csv:
        csv_fp = open(args.csv, 'wt')
        csv_w = csv.writer(csv_fp)

    for queryfile, query, query_moltype, query_ksize in loader:
        notify('loaded query: {}... (k={}, {})', query.name()[:30],
               query_ksize, query_moltype)

        results = []
        search_fn = SearchMinHashesFindBest().search

        # note, "ignore self" here may prevent using newer 'tree.search' fn.
        for leaf in tree.find(search_fn, query, args.threshold):
            if leaf.data.md5sum() != query.md5sum(): # ignore self.
                similarity = query.similarity(
                    leaf.data, ignore_abundance=args.ignore_abundance)
                results.append((similarity, leaf.data))

        best_hit_sim = 0.0
        best_hit_query_name = ""
        if results:
            results.sort(key=lambda x: -x[0])   # reverse sort on similarity
            best_hit_sim, best_hit_query = results[0]
            notify('for {}, found: {:.2f} {}', query.name(),
                                               best_hit_sim,
                                               best_hit_query.name())
            best_hit_query_name = best_hit_query.name()
        else:
            notify('for {}, no match found', query.name())

        if csv_w:
            csv_w.writerow([queryfile, query.name(), best_hit_query_name,
                           best_hit_sim])

    if loader.skipped_ignore:
        notify('skipped/ignore: {}', loader.skipped_ignore)
    if loader.skipped_nosig:
        notify('skipped/nosig: {}', loader.skipped_nosig)

    if csv_fp:
        csv_fp.close()


def gather(args):
    from .search import gather_databases, format_bp

    set_quiet(args.quiet, args.debug)
    moltype = sourmash_args.calculate_moltype(args)

    # load the query signature & figure out all the things
    query = sourmash_args.load_query_signature(args.query,
                                               ksize=args.ksize,
                                               select_moltype=moltype)
    notify('loaded query: {}... (k={}, {})', query.name()[:30],
                                             query.minhash.ksize,
                                             sourmash_args.get_moltype(query))

    # verify signature was computed right.
    if query.minhash.scaled == 0:
        error('query signature needs to be created with --scaled')
        sys.exit(-1)

    # downsample if requested
    if args.scaled:
        notify('downsampling query from scaled={} to {}',
               query.minhash.scaled, int(args.scaled))
        query.minhash = query.minhash.downsample_scaled(args.scaled)

    # empty?
    if not len(query.minhash):
        error('no query hashes!? exiting.')
        sys.exit(-1)

    # set up the search databases
    databases = sourmash_args.load_dbs_and_sigs(args.databases, query, False,
                                                 args.traverse_directory)

    if not len(databases):
        error('Nothing found to search!')
        sys.exit(-1)

    found = []
    weighted_missed = 1
    new_max_hash = query.minhash.max_hash
    next_query = query

    for result, weighted_missed, new_max_hash, next_query in gather_databases(query, databases, args.threshold_bp, args.ignore_abundance):
        if not len(found):                # first result? print header.
            if query.minhash.track_abundance and not args.ignore_abundance:
                print_results("")
                print_results("overlap     p_query p_match avg_abund")
                print_results("---------   ------- ------- ---------")
            else:
                print_results("")
                print_results("overlap     p_query p_match")
                print_results("---------   ------- -------")


        # print interim result & save in a list for later use
        pct_query = '{:.1f}%'.format(result.f_unique_weighted*100)
        pct_genome = '{:.1f}%'.format(result.f_match*100)
        average_abund ='{:.1f}'.format(result.average_abund)
        name = result.match._display_name(40)

        if query.minhash.track_abundance and not args.ignore_abundance:
            print_results('{:9}   {:>7} {:>7} {:>9}    {}',
                      format_bp(result.intersect_bp), pct_query, pct_genome,
                      average_abund, name)
        else:
            print_results('{:9}   {:>7} {:>7}    {}',
                      format_bp(result.intersect_bp), pct_query, pct_genome,
                      name)
        found.append(result)

        if args.num_results and len(found) >= args.num_results:
            break


    # basic reporting
    print_results('\nfound {} matches total;', len(found))
    if args.num_results and len(found) == args.num_results:
        print_results('(truncated gather because --num-results={})',
                      args.num_results)

    print_results('the recovered matches hit {:.1f}% of the query',
           (1 - weighted_missed) * 100)
    print_results('')

    if found and args.output:
        fieldnames = ['intersect_bp', 'f_orig_query', 'f_match',
                      'f_unique_to_query', 'f_unique_weighted',
                      'average_abund', 'median_abund', 'std_abund', 'name',
                      'filename', 'md5', 'f_match_orig']

        with FileOutput(args.output, 'wt') as fp:
            w = csv.DictWriter(fp, fieldnames=fieldnames)
            w.writeheader()
            for result in found:
                d = dict(result._asdict())
                del d['match']                 # actual signature not in CSV.
                w.writerow(d)

    if found and args.save_matches:
        notify('saving all matches to "{}"', args.save_matches)
        with FileOutput(args.save_matches, 'wt') as fp:
            sig.save_signatures([ r.match for r in found ], fp)

    if args.output_unassigned:
        if not len(query.minhash):
            notify('no unassigned hashes! not saving.')
        else:
            notify('saving unassigned hashes to "{}"', args.output_unassigned)

            with_abundance = next_query.minhash.track_abundance
            e = MinHash(ksize=query.minhash.ksize, n=0, max_hash=new_max_hash,
                        track_abundance=with_abundance)
            if with_abundance:
                abunds = next_query.minhash.get_mins(with_abundance=True)
                e.set_abundances(abunds)
            else:
                e.add_many(next_query.minhash.get_mins())

            with FileOutput(args.output_unassigned, 'wt') as fp:
                sig.save_signatures([ sig.SourmashSignature(e) ], fp)


def multigather(args):
    "Gather many signatures against multiple databases."
    from .search import gather_databases, format_bp

    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    if not args.db:
        error('Error! must specify at least one database with --db')
        sys.exit(-1)

    if not args.query:
        error('Error! must specify at least one query signature with --query')
        sys.exit(-1)

    # flatten --db and --query
    args.db = [item for sublist in args.db for item in sublist]
    args.query = [item for sublist in args.query for item in sublist]

    query = sourmash_args.load_query_signature(args.query[0],
                                               ksize=args.ksize,
                                               select_moltype=moltype)
    # set up the search databases
    databases = sourmash_args.load_dbs_and_sigs(args.db, query, False,
                                                args.traverse_directory)

    if not len(databases):
        error('Nothing found to search!')
        sys.exit(-1)

    # run gather on all the queries.
    for queryfile in args.query:
        # load the query signature & figure out all the things
        query = sourmash_args.load_query_signature(queryfile,
                                                   ksize=args.ksize,
                                                   select_moltype=moltype)
        notify('loaded query: {}... (k={}, {})', query.name()[:30],
                                                 query.minhash.ksize,
                                                 sourmash_args.get_moltype(query))

        # verify signature was computed right.
        if query.minhash.max_hash == 0:
            error('query signature needs to be created with --scaled; skipping')
            continue

        # downsample if requested
        if args.scaled:
            notify('downsampling query from scaled={} to {}',
                   query.minhash.scaled, int(args.scaled))
            query.minhash = query.minhash.downsample_scaled(args.scaled)

        # empty?
        if not len(query.minhash):
            error('no query hashes!? skipping to next..')
            continue

        found = []
        weighted_missed = 1
        for result, weighted_missed, new_max_hash, next_query in gather_databases(query, databases, args.threshold_bp, args.ignore_abundance):
            if not len(found):                # first result? print header.
                if query.minhash.track_abundance and not args.ignore_abundance:
                    print_results("")
                    print_results("overlap     p_query p_match avg_abund")
                    print_results("---------   ------- ------- ---------")
                else:
                    print_results("")
                    print_results("overlap     p_query p_match")
                    print_results("---------   ------- -------")


            # print interim result & save in a list for later use
            pct_query = '{:.1f}%'.format(result.f_unique_weighted*100)
            pct_genome = '{:.1f}%'.format(result.f_match*100)
            average_abund ='{:.1f}'.format(result.average_abund)
            name = result.match._display_name(40)

            if query.minhash.track_abundance and not args.ignore_abundance:
                print_results('{:9}   {:>7} {:>7} {:>9}    {}',
                          format_bp(result.intersect_bp), pct_query, pct_genome,
                          average_abund, name)
            else:
                print_results('{:9}   {:>7} {:>7}    {}',
                          format_bp(result.intersect_bp), pct_query, pct_genome,
                          name)
            found.append(result)


        # basic reporting
        print_results('\nfound {} matches total;', len(found))

        print_results('the recovered matches hit {:.1f}% of the query',
               (1 - weighted_missed) * 100)
        print_results('')

        if not found:
            notify('nothing found... skipping.')
            continue

        output_base = os.path.basename(queryfile)
        output_csv = output_base + '.csv'

        fieldnames = ['intersect_bp', 'f_orig_query', 'f_match',
                      'f_unique_to_query', 'f_unique_weighted',
                      'average_abund', 'median_abund', 'std_abund', 'name',
                      'filename', 'md5', 'f_match_orig']
        with open(output_csv, 'wt') as fp:
            w = csv.DictWriter(fp, fieldnames=fieldnames)
            w.writeheader()
            for result in found:
                d = dict(result._asdict())
                del d['match']                 # actual signature not in CSV.
                w.writerow(d)

        output_matches = output_base + '.matches.sig'
        with open(output_matches, 'wt') as fp:
            outname = output_matches
            notify('saving all matches to "{}"', outname)
            sig.save_signatures([ r.match for r in found ], fp)

        output_unassigned = output_base + '.unassigned.sig'
        with open(output_unassigned, 'wt') as fp:
            if not found:
                notify('nothing found - entire query signature unassigned.')
            elif not len(query.minhash):
                notify('no unassigned hashes! not saving.')
            else:
                notify('saving unassigned hashes to "{}"', output_unassigned)

                e = MinHash(ksize=query.minhash.ksize, n=0, max_hash=new_max_hash)
                e.add_many(next_query.minhash.get_mins())
                sig.save_signatures([ sig.SourmashSignature(e) ], fp)

        # fini, next query!


def watch(args):
    "Build a signature from raw FASTA/FASTQ coming in on stdin, search."
    set_quiet(args.quiet)

    if args.input_is_protein and args.dna:
        notify('WARNING: input is protein, turning off nucleotide hashing.')
        args.dna = False
        args.protein = True

    if args.dna and args.protein:
        notify('ERROR: cannot use "watch" with both nucleotide and protein.')

    if args.dna:
        moltype = 'DNA'
        is_protein = False
        dayhoff = False
        hp = False
    elif args.protein:
        moltype = 'protein'
        is_protein = True
        dayhoff = False
        hp = False
    elif args.dayhoff:
        moltype = 'dayhoff'
        is_protein = True
        dayhoff = True
        hp = False
    else:
        moltype = 'hp'
        is_protein = True
        dayhoff = False
        hp = True

    tree = load_sbt_index(args.sbt_name)

    # check ksize from the SBT we are loading
    ksize = args.ksize
    if ksize is None:
        leaf = next(iter(tree.leaves()))
        tree_mh = leaf.data.minhash
        ksize = tree_mh.ksize

    E = MinHash(ksize=ksize, n=args.num_hashes, is_protein=is_protein, dayhoff=dayhoff, hp=hp)

    notify('Computing signature for k={}, {} from stdin', ksize, moltype)

    def do_search():
        results = []
        streamsig = sig.SourmashSignature(E, filename='stdin', name=args.name)
        for similarity, match, _ in tree.search(streamsig,
                                                threshold=args.threshold,
                                                best_only=True,
                                                ignore_abundance=True,
                                                do_containment=False):
            results.append((similarity, match))

        return results

    notify('reading sequences from stdin')
    screed_iter = screed.open(args.inp_file)
    watermark = WATERMARK_SIZE

    # iterate over input records
    n = 0
    for n, record in enumerate(screed_iter):
        # at each watermark, print status & check cardinality
        if n >= watermark:
            notify('\r... read {} sequences', n, end='')
            watermark += WATERMARK_SIZE

            if do_search():
                break

        if args.input_is_protein:
            E.add_protein(record.sequence)
        else:
            E.add_sequence(record.sequence, False)

    results = do_search()
    if not results:
        notify('... read {} sequences, no matches found.', n)
    else:
        results.sort(key=lambda x: -x[0])   # take best
        similarity, found_sig = results[0]
        print_results('FOUND: {}, at {:.3f}', found_sig.name(),
               similarity)

    if args.output:
        notify('saving signature to {}', args.output)
        with FileOutput(args.output, 'wt') as fp:
            streamsig = sig.SourmashSignature(E, filename='stdin',
                                              name=args.name)
        sig.save_signatures([streamsig], fp)


def migrate(args):
    "Migrate an SBT database to the latest version."
    tree = load_sbt_index(args.sbt_name, print_version_warning=False)

    notify('saving SBT under "{}".', args.sbt_name)
    tree.save(args.sbt_name, structure_only=True)
