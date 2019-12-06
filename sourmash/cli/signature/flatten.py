from argparse import FileType
import sourmash
from sourmash.cli.utils import add_moltype_args, add_ksize_arg

def subparser(subparsers):
    subparser = subparsers.add_parser('flatten')
    subparser.add_argument('signatures', nargs='+')
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    subparser.add_argument(
        '-o', '--output', metavar='FILE', type=FileType('wt'),
        help='output signature to this file'
    )
    subparser.add_argument(
        '--md5', default=None,
        help='select signatures whose md5 contains this substring'
    )
    subparser.add_argument(
        '--name', default=None,
        help='select signatures whose name contains this substring'
    )
    add_ksize_arg(subparser, 31)
    add_moltype_args(subparser)
