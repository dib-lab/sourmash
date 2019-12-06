from argparse import FileType
import sourmash
from sourmash.cli.utils import add_moltype_args, add_ksize_arg

def subparser(subparsers):
    subparser = subparsers.add_parser('downsample')
    subparser.add_argument('signatures', nargs="+")
    subparser.add_argument(
        '--scaled', type=int, default=0,
        help='scaled value to downsample to'
    )
    subparser.add_argument(
        '--num', metavar='N', type=int, default=0,
        help='num value to downsample to'
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE', type=FileType('wt'),
        help='output signature to this file'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)
