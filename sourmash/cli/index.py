"""index signatures for rapid search"""

usage="""

   sourmash index -k 31 dbname *.sig

Create an on-disk database of signatures that can be searched in low
memory with 'search' and 'gather'. All signatures must be the same
k-mer size, molecule type, and num/scaled; the standard signature
selectors (-k/--ksize, --scaled, --dna/--protein) choose which
signatures to be added.

The key options for index are:

 * `-k/--ksize <int>`: k-mer size to select
 * `--dna` or --protein`: nucleotide or protein signatures (default `--dna`)

If `dbname` ends with `.sbt.json`, index will create the database as a
collection of multiple files, with an index `dbname.sbt.json` and a
subdirectory `.sbt.dbname`. If `dbname` ends with `.sbt.zip`, index
will create a zip archive containing the multiple files. For sourmash
v2 and v3, `sbt.json` will be added automatically; this behavior will
change in sourmash v4 to default to `.sbt.zip`.

---
"""

from sourmash.cli.utils import add_moltype_args, add_ksize_arg


def subparser(subparsers):
    subparser = subparsers.add_parser('index', description=__doc__,
                                      usage=usage)
    subparser.add_argument('sbt_name', help='name to save index into; .sbt.zip or .sbt.json file')
    subparser.add_argument(
        'signatures', nargs='*',
        help='signatures to load into SBT'
    )
    subparser.add_argument(
        '--from-file',
        help='a file containing a list of signatures file to load'
    )
    subparser.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    add_ksize_arg(subparser, 31)
    subparser.add_argument(
        '-d', '--n_children', metavar='D', type=int, default=2,
        help='number of children for internal nodes; default=2'
    )
    subparser.add_argument(
        '--append', action='store_true', default=False,
        help='add signatures to an existing SBT'
    )
    subparser.add_argument(
        '-x', '--bf-size', metavar='S', type=float, default=1e5,
        help='Maximum Bloom filter size used for internal nodes. Set to 0 '
             'to calculate an optimal size given the complexity of the datasets '
             '(using a HyperLogLog to estimate the number of unique k-mers'
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='try loading *all* files in provided subdirectories, not just .sig files"'
    )
    subparser.add_argument(
        '-s', '--sparseness', metavar='FLOAT', type=float, default=.0,
        help='What percentage of internal nodes will not be saved; ranges '
        'from 0.0 (save all nodes) to 1.0 (no nodes saved)'
    )
    subparser.add_argument(
        '--scaled', metavar='FLOAT', type=float, default=0,
        help='downsample signatures to the specified scaled factor'
    )
    add_moltype_args(subparser)


def main(args):
    import sourmash
    return sourmash.commands.index(args)
