"""search a metagenome signature against dbs"""

from sourmash.cli.utils import add_ksize_arg, add_moltype_args


def subparser(subparsers):
    subparser = subparsers.add_parser('gather')
    subparser.add_argument('query', help='query signature')
    subparser.add_argument(
        'databases', nargs='+',
        help='signatures/SBTs to search',
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true'
    )
    subparser.add_argument(
        '-n', '--num-results', default=None, type=int, metavar='N',
        help='number of results to report (default: terminate at --threshold-bp)'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output CSV containing matches to this file'
    )
    subparser.add_argument(
        '--save-matches', metavar='FILE',
        help='save gather matched signatures from the database to the '
        'specified file'
    )
    subparser.add_argument(
        '--save-prefetch', metavar='FILE',
        help='save all prefetch-matched signatures from the databases to the '
        'specified file or directory'
    )
    subparser.add_argument(
        '--threshold-bp', metavar='REAL', type=float, default=5e4,
        help='reporting threshold (in bp) for estimated overlap with remaining query (default=50kb)'
    )
    subparser.add_argument(
        '--output-unassigned', metavar='FILE',
        help='output unassigned portions of the query as a signature to the '
        'specified file'
    )
    subparser.add_argument(
        '--scaled', metavar='FLOAT', type=float, default=0,
        help='downsample query to the specified scaled factor'
    )
    subparser.add_argument(
        '--ignore-abundance',  action='store_true',
        help='do NOT use k-mer abundances if present'
    )
    subparser.add_argument(
        '--md5', default=None,
        help='select the signature with this md5 as query'
    )
    subparser.add_argument(
        '--cache-size', default=0, type=int, metavar='N',
        help='number of internal SBT nodes to cache in memory (default: 0, cache all nodes)'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)

    # advanced parameters
    subparser.add_argument(
        '--linear', dest="linear", action='store_true',
        help="force a low-memory but maybe slower database search",
    )
    subparser.add_argument(
        '--no-linear', dest="linear", action='store_false',
    )
    subparser.add_argument(
        '--no-prefetch', dest="prefetch", action='store_false',
        help="do not use prefetch before gather; see documentation",
    )
    subparser.add_argument(
        '--prefetch', dest="prefetch", action='store_true',
        help="use prefetch before gather; see documentation",
    )


def main(args):
    import sourmash
    return sourmash.commands.gather(args)
