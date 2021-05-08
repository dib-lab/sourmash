"""search a signature against dbs, find all overlaps"""

from sourmash.cli.utils import add_ksize_arg, add_moltype_args


def subparser(subparsers):
    subparser = subparsers.add_parser('prefetch')
    subparser.add_argument('query', help='query signature')
    subparser.add_argument("databases",
        nargs="*",
        help="one or more databases to search",
    )
    subparser.add_argument(
        "--db-from-file",
        default=None,
        help="list of paths containing signatures to search"
    )
    subparser.add_argument(
        "--linear", action='store_true',
        help="force linear traversal of indexes to minimize loading time and memory use"
    )
    subparser.add_argument(
        '--no-linear', dest="linear", action='store_false',
    )

    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-d', '--debug', action='store_true'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output CSV containing matches to this file'
    )
    subparser.add_argument(
        '--save-matches', metavar='FILE',
        help='save all matching signatures from the databases to the '
        'specified file or directory'
    )
    subparser.add_argument(
        '--threshold-bp', metavar='REAL', type=float, default=5e4,
        help='reporting threshold (in bp) for estimated overlap with remaining query hashes (default=50kb)'
    )
    subparser.add_argument(
        '--save-unmatched-hashes', metavar='FILE',
        help='output unmatched query hashes as a signature to the '
        'specified file'
    )
    subparser.add_argument(
        '--save-matching-hashes', metavar='FILE',
        help='output matching query hashes as a signature to the '
        'specified file'
    )
    subparser.add_argument(
        '--scaled', metavar='FLOAT', type=float, default=None,
        help='downsample signatures to the specified scaled factor'
    )
    subparser.add_argument(
        '--md5', default=None,
        help='select the signature with this md5 as query'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.commands.prefetch(args)
