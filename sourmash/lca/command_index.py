#! /usr/bin/env python
"""
Build a lowest-common-ancestor database with given taxonomy and genome sigs.
"""
import sys
import csv
from collections import defaultdict
import os.path

from sourmash import sourmash_args, load_signatures
from sourmash.sourmash_args import load_file_as_signatures
from sourmash.logging import notify, error, debug, set_quiet
from . import lca_utils
from .lca_utils import LineagePair
from .lca_db import LCA_Database
from sourmash.sourmash_args import DEFAULT_LOAD_K


def load_taxonomy_assignments(filename, delimiter=',', start_column=2,
                              use_headers=True, force=False):
    """
    Load a taxonomy assignment spreadsheet into a dictionary.

    The 'assignments' dictionary that's returned maps identifiers to
    lineage tuples.
    """
    mode = 'rt'

    # parse spreadsheet!
    fp = open(filename, mode)
    r = csv.reader(fp, delimiter=delimiter)
    row_headers = ['identifiers']
    row_headers += ['_skip_']*(start_column - 2)
    row_headers += list(lca_utils.taxlist())

    # first check that headers are interpretable.
    if use_headers:
        notify('examining spreadsheet headers...')
        first_row = next(iter(r))

        n_disagree = 0
        for (column, value) in zip(row_headers, first_row):
            if column == '_skip_':
                continue

            if column.lower() != value.lower():
                notify("** assuming column '{}' is {} in spreadsheet",
                       value, column)
                n_disagree += 1
                if n_disagree > 2:
                    error('whoa, too many assumptions. are the headers right?')
                    error('expecting {}', ",".join(row_headers))
                    if not force:
                        sys.exit(-1)
                    notify('...continue, because --force was specified.')

    # convert into a lineage pair
    assignments = {}
    num_rows = 0
    n_species = 0
    n_strains = 0
    for row in r:
        if row and row[0].strip():        # want non-empty row
            num_rows += 1
            lineage = list(zip(row_headers, row))
            lineage = [ x for x in lineage if x[0] != '_skip_' ]

            ident = lineage[0][1]
            lineage = lineage[1:]

            # clean lineage of null names, replace with 'unassigned'
            lineage = [ (a, lca_utils.filter_null(b)) for (a,b) in lineage ]
            lineage = [ LineagePair(a, b) for (a, b) in lineage ]

            # remove end nulls
            while lineage and lineage[-1].name == 'unassigned':
                lineage = lineage[:-1]

            # store lineage tuple
            if lineage:
                # check duplicates
                if ident in assignments:
                    if assignments[ident] != tuple(lineage):
                        if not force:
                            raise Exception("multiple lineages for identifier {}".format(ident))
                else:
                    assignments[ident] = tuple(lineage)

                    if lineage[-1].rank == 'species':
                        n_species += 1
                    elif lineage[-1].rank == 'strain':
                        n_species += 1
                        n_strains += 1

    fp.close()

    # this is to guard against a bug that happened once and I can't find
    # any more, when building a large GTDB-based database :) --CTB
    if len(assignments) * 0.2 > n_species and len(assignments) > 50:
        if not force:
            error('')
            error("ERROR: fewer than 20% of lineages have species-level resolution!?")
            error("({} species assignments found, of {} assignments total)",
                  n_species, len(assignments))
            error("** If this is intentional, re-run the command with -f.")
            sys.exit(-1)

    return assignments, num_rows


def generate_report(record_duplicates, record_no_lineage, record_remnants,
                    unused_lineages, unused_identifiers, filename):
    """
    Output a report of anomalies from building the index.
    """
    with open(filename, 'wt') as fp:
        print('Duplicate signatures:', file=fp)
        fp.write("\n".join(record_duplicates))
        fp.write("\n")
        print('----\nUnused identifiers:', file=fp)
        fp.write("\n".join(unused_identifiers))
        fp.write("\n")
        print('----\nNo lineage provided for these identifiers:', file=fp)
        fp.write("\n".join(record_no_lineage))
        fp.write("\n")
        print('----\nNo signatures found for these identifiers:', file=fp)
        fp.write('\n'.join(record_remnants))
        fp.write("\n")
        print('----\nUnused lineages:', file=fp)
        for lineage in unused_lineages:
            fp.write(";".join(lca_utils.zip_lineage(lineage)))
            fp.write("\n")


def index(args):
    """
    main function for building an LCA database.
    """
    if args.start_column < 2:
        error('error, --start-column cannot be less than 2')
        sys.exit(-1)

    set_quiet(args.quiet, args.debug)

    args.scaled = int(args.scaled)

    if args.ksize is None:
        args.ksize = DEFAULT_LOAD_K

    moltype = sourmash_args.calculate_moltype(args, default='DNA')

    notify('Building LCA database with ksize={} scaled={} moltype={}.',
           args.ksize, args.scaled, moltype)

    # first, load taxonomy spreadsheet
    delimiter = ','
    if args.tabs:
        delimiter = '\t'
    assignments, num_rows = load_taxonomy_assignments(args.csv,
                                               delimiter=delimiter,
                                               start_column=args.start_column,
                                               use_headers=not args.no_headers,
                                               force=args.force)

    notify('{} distinct identities in spreadsheet out of {} rows.',
           len(assignments), num_rows)
    notify('{} distinct lineages in spreadsheet out of {} rows.',
           len(set(assignments.values())), num_rows)

    db = LCA_Database(args.ksize, args.scaled, moltype)

    inp_files = list(args.signatures)
    if args.from_file:
        more_files = sourmash_args.load_file_list_of_signatures(args.from_file)
        inp_files.extend(more_files)

    # track duplicates
    md5_to_name = {}

    #
    # main loop, connecting lineage ID to signature.
    #

    n = 0
    total_n = len(inp_files)
    record_duplicates = set()
    record_no_lineage = set()
    record_remnants = set(assignments)
    record_used_lineages = set()
    record_used_idents = set()
    n_skipped = 0
    for filename in inp_files:
        n += 1
        it = load_file_as_signatures(filename, ksize=args.ksize,
                                     select_moltype=moltype,
                                     yield_all_files=args.force)
        for sig in it:
            notify(u'\r\033[K', end=u'')
            notify('\r... loading signature {} ({} of {}); skipped {} so far', str(sig)[:30], n, total_n, n_skipped, end='')
            debug(filename, sig)

            # block off duplicates.
            if sig.md5sum() in md5_to_name:
                debug('WARNING: in file {}, duplicate md5sum: {}; skipping', filename, sig.md5sum())
                record_duplicates.add(filename)
                continue

            md5_to_name[sig.md5sum()] = str(sig)

            # parse identifier, potentially with splitting
            ident = sig.name
            if args.split_identifiers: # hack for NCBI-style names, etc.
                # split on space...
                ident = ident.split(' ')[0]
                # ...and on period.
                ident = ident.split('.')[0]
            elif args.use_filename_as_identifier:
                ident = os.path.basename(sig.filename)

            lineage = assignments.get(ident)

            # punt if no lineage and --require-taxonomy
            if lineage is None and args.require_taxonomy:
                debug('(skipping, because --require-taxonomy was specified)')
                n_skipped += 1
                continue

            # add the signature into the database.
            try:
                db.insert(sig, ident=ident, lineage=lineage)
            except ValueError as e:
                error("ERROR: cannot insert signature '{}' (md5 {}, loaded from '{}') into database.",
                      sig, sig.md5sum()[:8], filename)
                error("ERROR: {}", str(e))
                sys.exit(-1)

            if lineage:
                # remove from our list of remaining ident -> lineage
                record_remnants.remove(ident)

                # track ident as used
                record_used_idents.add(ident)
                record_used_lineages.add(lineage)

            # track lineage info - either no lineage, or this lineage used.
            else:
                debug('WARNING: no lineage assignment for {}.', ident)
                record_no_lineage.add(ident)

    # end main add signatures loop

    if n_skipped:
        notify('... loaded {} signatures; skipped {} because of --require-taxonomy.', total_n, n_skipped)
    else:
        notify('... loaded {} signatures.', total_n)

    # check -- did we find any signatures?
    if n == 0:
        error('ERROR: no signatures found. ??')
        sys.exit(1)

    # check -- did the signatures we found have any hashes?
    if not db.hashval_to_idx:
        error('ERROR: no hash values found - are there any signatures?')
        sys.exit(1)
    notify('loaded {} hashes at ksize={} scaled={}', len(db.hashval_to_idx),
           args.ksize, args.scaled)

    # summarize:
    notify('{} assigned lineages out of {} distinct lineages in spreadsheet.',
           len(record_used_lineages), len(set(assignments.values())))
    unused_lineages = set(assignments.values()) - record_used_lineages

    notify('{} identifiers used out of {} distinct identifiers in spreadsheet.',
           len(record_used_idents), len(set(assignments)))

    assert record_used_idents.issubset(set(assignments))
    unused_identifiers = set(assignments) - record_used_idents

    # now, save!
    db_outfile = args.lca_db_out
    if not (db_outfile.endswith('.lca.json') or \
                db_outfile.endswith('.lca.json.gz')):   # logic -> db.save
        db_outfile += '.lca.json'
    notify('saving to LCA DB: {}'.format(db_outfile))

    db.save(db_outfile)

    ## done!

    # output a record of stuff if requested/available:
    if record_duplicates or record_no_lineage or record_remnants or unused_lineages:
        if record_duplicates:
            notify('WARNING: {} duplicate signatures.', len(record_duplicates))
        if record_no_lineage:
            notify('WARNING: no lineage provided for {} signatures.',
                   len(record_no_lineage))
        if record_remnants:
            notify('WARNING: no signatures for {} spreadsheet rows.',
                   len(record_remnants))
        if unused_lineages:
            notify('WARNING: {} unused lineages.', len(unused_lineages))

        if unused_identifiers:
            notify('WARNING: {} unused identifiers.', len(unused_identifiers))

        if args.report:
            notify("generating a report and saving in '{}'", args.report)
            generate_report(record_duplicates, record_no_lineage,
                            record_remnants, unused_lineages,
                            unused_identifiers, args.report)
        else:
            notify('(You can use --report to generate a detailed report.)')


if __name__ == '__main__':
    sys.exit(index(sys.argv[1:]))
