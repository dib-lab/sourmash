"""Define the command line interface for sourmash sig

The top level CLI is defined in ../__init__.py. This module defines the CLI for
`sourmash sig` operations.
"""

from . import cat
from . import split
from . import describe
from . import downsample
from . import extract
from . import filter
from . import flatten
from . import abundhist
from . import intersect
from . import merge
from . import rename
from . import subtract
from . import ingest
from . import export
from . import overlap
from ..utils import command_list
from argparse import SUPPRESS, RawDescriptionHelpFormatter
import os
import sys


def subparser(subparsers):
    subparser = subparsers.add_parser('sig', formatter_class=RawDescriptionHelpFormatter, usage=SUPPRESS, aliases=['signature'])
    desc = 'Operations\n'
    clidir = os.path.dirname(__file__)
    ops = command_list(clidir)
    for subcmd in ops:
        docstring = getattr(sys.modules[__name__], subcmd).__doc__
        helpstring = 'sourmash sig {op:s} --help'.format(op=subcmd)
        desc += '        {hs:33s} {ds:s}\n'.format(hs=helpstring, ds=docstring)
    s = subparser.add_subparsers(
        title='Manipulate signature files', dest='subcmd', metavar='subcmd', help=SUPPRESS,
        description=desc
    )
    for subcmd in ops:
        getattr(sys.modules[__name__], subcmd).subparser(s)
    subparser._action_groups.reverse()
    subparser._optionals.title = 'Options'
