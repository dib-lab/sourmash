from argparse import ArgumentParser, RawDescriptionHelpFormatter, SUPPRESS
import os
import sys

import sourmash

from . import utils

# Commands
from . import categorize
from . import compare_csv
from . import compute
from . import dump
from . import gather
from . import import_csv
from . import info
from . import index
from . import migrate
from . import multigather
from . import plot
from . import sbt_combine
from . import search
from . import watch

# Subcommand groups
from . import lca
from . import sig
from . import storage


class SourmashParser(ArgumentParser):
    def __init__(self, citation=True, **kwargs):
        super(SourmashParser, self).__init__(**kwargs)
        self.citation = citation
        self._citation_printed = False

    def print_citation(self):
        if self._citation_printed:
            return
        from sourmash.logging import notify
        notify("== This is sourmash version {version}. ==", version=sourmash.VERSION)
        notify("== Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==\n")
        self._citation_printed = True

    def _subparser_from_name(self, name):
        """Given a name, get the subparser instance registered with this parser."""
        container = self._actions
        if name is None:
            return None
        for action in container:
            if action.choices is None:
                continue
            elif name in action.choices:
                return action.choices[name]

    def parse_args(self, args=None, namespace=None):
        if (args is None and len(sys.argv) == 1) or (args is not None and len(args) == 0):
            self.print_help()
            raise SystemExit(1)
        args = super(SourmashParser, self).parse_args(args=args, namespace=namespace)
        if ('quiet' not in args or not args.quiet) and self.citation:
            self.print_citation()

        if 'subcmd' in args and args.subcmd is None:
            self._subparser_from_name(args.cmd).print_help()
            raise SystemExit(1)

        # BEGIN: dirty hacks to simultaneously support new and previous interface
        if hasattr(args, 'subcmd') and args.subcmd == 'import':
            args.subcmd = 'ingest'
        # END: dirty hacks to simultaneously support new and previous interface
        return args


def get_parser():
    module_descs = {
        'lca': 'Lowest common ancestor (LCA) based operations',
        'sig': 'Operations on signatures (MinHash sketches)',
        'storage': 'Operations on storage',
    }

    clidir = os.path.dirname(__file__)
    basic_ops = utils.command_list(clidir)
    usage = '    Basic operations\n'
    for bo in basic_ops:
        usage += '        sourmash {op:s} --help\n'.format(op=bo)
    cmd_group_dirs = next(os.walk(clidir))[1]
    cmd_group_dirs = filter(utils.opfilter, cmd_group_dirs)
    cmd_group_dirs = sorted(cmd_group_dirs)
    for dirpath in cmd_group_dirs:
        usage += '\n    ' + module_descs[dirpath] + '\n'
        usage += '        sourmash {gd:s} --help\n'.format(gd=dirpath)

    desc = 'Compute, compare, manipulate, and analyze MinHash sketches of DNA sequences.\n\nUsage instructions:\n' + usage
    parser = SourmashParser(prog='sourmash', description=desc, formatter_class=RawDescriptionHelpFormatter, usage=SUPPRESS)
    parser._optionals.title = 'Options'
    parser.add_argument('-v', '--version', action='version', version='sourmash '+ sourmash.VERSION)
    parser.add_argument('-q', '--quiet', action='store_true', help='don\'t print citation information')
    sub = parser.add_subparsers(
        title='Instructions', dest='cmd', metavar='cmd', help=SUPPRESS,
    )
    for op in basic_ops + cmd_group_dirs:
        getattr(sys.modules[__name__], op).subparser(sub)
    parser._action_groups.reverse()
    return parser
