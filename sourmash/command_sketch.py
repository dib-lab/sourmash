"""
Functions implementing the 'sketch' subcommands and related functions.
"""
import sys

from .signature import SourmashSignature
from .logging import notify, error, set_quiet
from .command_compute import (_compute_individual, _compute_merged,
                              ComputeParameters)


DEFAULTS = dict(
    dna='k=31,scaled=1000,noabund',
    protein='k=21,scaled=200,noabund',
    dayhoff='k=19,scaled=200,noabund',
    hp='k=30,scaled=200,noabund'
)


def _parse_params_str(params_str):
    "Parse a parameter string of the form 'k=ks,num=num,scaled=scaled,abund'."
    moltype = None
    d = {}
    pp = params_str.split(',')
    for p in pp:
        if p == 'abund':
            d['track_abundance'] = True
        elif p == 'noabund':
            d['track_abundance'] = False
        elif p.startswith('k='):
            d['ksize'] = int(p[2:])
        elif p.startswith('num='):
            if d.get('scaled'):
                raise ValueError("cannot set both num and scaled in a single minhash")
            d['num'] = int(p[4:])
            d['scaled'] = 0
        elif p.startswith('scaled='):
            if d.get('num'):
                raise ValueError("cannot set both num and scaled in a single minhash")
            d['scaled'] = int(p[7:])
            d['num'] = 0
        elif p.startswith('seed='):
            d['seed'] = int(p[5:])
        elif p == 'protein':
            moltype = 'protein'
        elif p == 'dayhoff':
            moltype = 'dayhoff'
        elif p == 'hp':
            moltype = 'hp'
        elif p == 'dna':
            moltype = 'dna'
        else:
            raise ValueError(f"unknown component '{p}' in params string")

    return moltype, d


class _signatures_for_sketch_factory(object):
    "Build sigs on demand, based on args input to 'sketch'."
    def __init__(self, params_str_list, default_moltype, mult_ksize_by_3):
        defaults = {}
        for moltype, pstr in DEFAULTS.items():
            mt, d = _parse_params_str(pstr)
            assert mt is None
            defaults[moltype] = d
        self.defaults = defaults

        self.params_list = []
        self.mult_ksize_by_3 = mult_ksize_by_3
        for params_str in params_str_list:
            moltype, d = _parse_params_str(params_str)
            if moltype is None:
                moltype = default_moltype
            self.params_list.append((moltype, d))

    def __call__(self):
        "Produce a new set of signatures built to match the param strings."
        x = []

        for moltype, d in self.params_list:
            z = self.defaults[moltype]
            def_ksize = z['ksize']
            def_seed = z.get('seed', 42)
            def_num = z.get('num', 0)
            def_abund = z['track_abundance']
            def_scaled = z.get('scaled', 0)
            def_dna = z.get('is_dna', moltype == 'dna')
            def_protein = z.get('is_protein',  moltype == 'protein')
            def_dayhoff = z.get('is_dayhoff', moltype == 'dayhoff')
            def_hp = z.get('is_hp', moltype == 'hp')

            ksize = int(d.get('ksize', def_ksize))
            if self.mult_ksize_by_3:
                ksize = ksize*3

            params = ComputeParameters([ksize],
                                       d.get('seed', def_seed),
                                       def_protein,
                                       def_dayhoff,
                                       def_hp,
                                       def_dna,
                                       d.get('num', def_num),
                                       d.get('track_abundance', def_abund),
                                       d.get('scaled', def_scaled))
            sig = SourmashSignature.from_params(params)
            x.append(sig)

        return x


def _execute_sketch(args, signatures_factory):
    "Once configured, run 'sketch' the same way underneath."
    set_quiet(args.quiet)

    if args.license != 'CC0':
        error('error: sourmash only supports CC0-licensed signatures. sorry!')
        sys.exit(-1)

    notify('computing signatures for files: {}', ", ".join(args.filenames))

    if args.merge and not args.output:
        error("ERROR: must specify -o with --merge")
        sys.exit(-1)

    if args.output and args.outdir:
        error("ERROR: --outdir doesn't make sense with -o/--output")
        sys.exit(-1)

    # get number of output sigs:
    num_sigs = len(signatures_factory.params_list)
    notify('Computing a total of {} signature(s).', num_sigs)

    if num_sigs == 0:
        error('...nothing to calculate!? Exiting!')
        sys.exit(-1)

    if args.merge:               # single name specified - combine all
        _compute_merged(args, signatures_factory)
    else:                        # compute individual signatures
        _compute_individual(args, signatures_factory)


def dna(args):
    """Compute a DNA signature for one or more files.

    CTB: make usable via Python?
    """
    # TBD/FIXME
    args.input_is_10x = False # CTB

    # for dna:
    args.input_is_protein = False

    # provide good defaults for dna
    if not args.param_string:
        args.param_string = ['k=31,scaled=1000,noabund']

    signatures_factory = _signatures_for_sketch_factory(args.param_string,
                                                        'dna',
                                                        mult_ksize_by_3=False)

    _execute_sketch(args, signatures_factory)


def protein(args):
    """Compute a protein signature for one or more files.

    CTB: make usable via Python?
    """
    # for protein:
    args.input_is_10x = False
    args.input_is_protein = True

    # provide good defaults for dayhoff/hp/protein!
    if args.dayhoff and args.hp:
        raise ValueError("cannot set both --dayhoff and --hp")
    if args.dayhoff:
        moltype = 'dayhoff'
        default_param_string = ['k=19,scaled=200,noabund']
    elif args.hp:
        moltype = 'hp'
        default_param_string = ['k=30,scaled=200,noabund']
    else:
        moltype = 'protein'
        default_param_string = ['k=21,scaled=200,noabund']

    if not args.param_string:
        args.param_string = default_param_string

    signatures_factory = _signatures_for_sketch_factory(args.param_string,
                                                        moltype,
                                                        mult_ksize_by_3=True)

    _execute_sketch(args, signatures_factory)


def translate(args):
    """Compute protein signatures from DNA/RNA, for one or more files.

    CTB: make usable via Python?
    """
    # for translate:
    args.input_is_10x = False
    args.input_is_protein = False

    # provide good defaults for dayhoff/hp/protein!
    if args.dayhoff and args.hp:
        raise ValueError("cannot set both --dayhoff and --hp")
    if args.dayhoff:
        moltype = 'dayhoff'
        default_param_string = ['k=19,scaled=200,noabund']
    elif args.hp:
        moltype = 'hp'
        default_param_string = ['k=30,scaled=200,noabund']
    else:
        moltype = 'protein'
        default_param_string = ['k=21,scaled=200,noabund']

    if not args.param_string:
        args.param_string = default_param_string

    signatures_factory = _signatures_for_sketch_factory(args.param_string,
                                                        moltype,
                                                        mult_ksize_by_3=True)

    _execute_sketch(args, signatures_factory)
