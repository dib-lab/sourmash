"""output abundance histogram"""

import sys

from sourmash.cli.utils import add_moltype_args, add_ksize_arg


def subparser(subparsers):
    subparser = subparsers.add_parser('abundhist')
    subparser.add_argument('signatures', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE',
        help='output histogram to this file (in CSV format)'
    )
    subparser.add_argument(
        '--abundances', metavar='FILE',
        help='output hashes and abundances to this file (in CSV format)')
    subparser.add_argument(
        '--md5', default=None,
        help='select signatures whose md5 contains this substring'
    )
    subparser.add_argument(
        '--name', default=None,
        help='select signatures whose name contains this substring'
    )
    subparser.add_argument(
        '--max', type=int, default=None,
        help='max value for histogram range (default none)')
    subparser.add_argument(
        '--min', type=int, default=None,
        help='min value for histogram range (default none)')
    subparser.add_argument(
        '--bins', type=int, default=10,
        help='number of bins (default 10)')
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.sig.__main__.abundhist(args)
