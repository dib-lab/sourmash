"""
Functions implementing the 'compute' command and related functions.
"""
import os
import os.path
import sys
import random
import screed
import time

from . import sourmash_args
from .signature import SourmashSignature, save_signatures
from .logging import notify, error, set_quiet, short_notify
from .utils import RustObject
from ._lowlevel import ffi, lib

DEFAULT_COMPUTE_K = '21,31,51'
DEFAULT_LINE_COUNT = 1500


def compute(args):
    """Compute the signature for one or more files.

    Use cases:
        sourmash compute multiseq.fa              => multiseq.fa.sig, etc.
        sourmash compute genome.fa --singleton    => genome.fa.sig
        sourmash compute file1.fa file2.fa -o file.sig
            => creates one output file file.sig, with one signature for each
               input file.
        sourmash compute file1.fa file2.fa --merge merged -o file.sig
            => creates one output file file.sig, with all sequences from
               file1.fa and file2.fa combined into one signature.
    """
    set_quiet(args.quiet)

    if args.license != 'CC0':
        error('error: sourmash only supports CC0-licensed signatures. sorry!')
        sys.exit(-1)

    if args.input_is_protein and args.dna:
        notify('WARNING: input is protein, turning off nucleotide hashing')
        args.dna = False
        args.protein = True

    if args.scaled:
        if args.scaled < 1:
            error('ERROR: --scaled value must be >= 1')
            sys.exit(-1)
        if args.scaled != round(args.scaled, 0):
            error('ERROR: --scaled value must be integer value')
            sys.exit(-1)
        if args.scaled >= 1e9:
            notify('WARNING: scaled value is nonsensical!? Continuing anyway.')

        if args.num_hashes != 0:
            notify('setting num_hashes to 0 because --scaled is set')
            args.num_hashes = 0

    notify('computing signatures for files: {}', ", ".join(args.filenames))

    if args.randomize:
        notify('randomizing file list because of --randomize')
        random.shuffle(args.filenames)

    # get list of k-mer sizes for which to compute sketches
    ksizes = args.ksizes

    notify('Computing signature for ksizes: {}', str(ksizes))
    num_sigs = 0
    if args.dna and args.protein:
        notify('Computing both nucleotide and protein signatures.')
        num_sigs = 2*len(ksizes)
    elif args.dna and args.dayhoff:
        notify('Computing both nucleotide and Dayhoff-encoded protein '
               'signatures.')
        num_sigs = 2*len(ksizes)
    elif args.dna and args.hp:
        notify('Computing both nucleotide and hp-encoded protein '
               'signatures.')
        num_sigs = 2*len(ksizes)
    elif args.dna:
        notify('Computing only nucleotide (and not protein) signatures.')
        num_sigs = len(ksizes)
    elif args.protein:
        notify('Computing only protein (and not nucleotide) signatures.')
        num_sigs = len(ksizes)
    elif args.dayhoff:
        notify('Computing only Dayhoff-encoded protein (and not nucleotide) '
               'signatures.')
        num_sigs = len(ksizes)
    elif args.hp:
        notify('Computing only hp-encoded protein (and not nucleotide) '
               'signatures.')
        num_sigs = len(ksizes)

    if (args.protein or args.dayhoff or args.hp):
        notify("")
        notify("WARNING: you are using 'compute' to make a protein/dayhoff/hp signature,")
        notify("WARNING: but the meaning of ksize has changed in 4.0. Please see the")
        notify("WARNING: migration guide to sourmash v4.0 at http://sourmash.rtfd.io/")
        notify("")
        bad_ksizes = [ str(k) for k in ksizes if k % 3 != 0 ]
        if bad_ksizes:
            error('protein ksizes must be divisible by 3, sorry!')
            error('bad ksizes: {}', ", ".join(bad_ksizes))
            sys.exit(-1)

    notify('Computing a total of {} signature(s).', num_sigs)

    if num_sigs == 0:
        error('...nothing to calculate!? Exiting!')
        sys.exit(-1)

    if args.merge and not args.output:
        error("ERROR: must specify -o with --merge")
        sys.exit(-1)

    if args.output and args.outdir:
        error("ERROR: --outdir doesn't make sense with -o/--output")
        sys.exit(-1)

    if args.track_abundance:
        notify('Tracking abundance of input k-mers.')

    signatures_factory = _signatures_for_compute_factory(args)

    if args.merge:               # single name specified - combine all
        _compute_merged(args, signatures_factory)
    else:                        # compute individual signatures
        _compute_individual(args, signatures_factory)


class _signatures_for_compute_factory(object):
    "Build signatures on demand, based on args input to 'compute'."
    def __init__(self, args):
        self.args = args

    def __call__(self):
        args = self.args
        params = ComputeParameters(args.ksizes, args.seed, args.protein,
                                   args.dayhoff, args.hp, args.dna,
                                   args.num_hashes,
                                   args.track_abundance, args.scaled)
        sig = SourmashSignature.from_params(params)
        return [sig]


def _compute_individual(args, signatures_factory):
    siglist = []

    for filename in args.filenames:
        sigfile = os.path.basename(filename) + '.sig'
        if args.outdir:
            sigfile = os.path.join(args.outdir, sigfile)

        if not args.output and os.path.exists(sigfile) and not \
            args.force:
            short_notify('skipping {} - already done', filename)
            continue

        if args.singleton:
            siglist = []
            for n, record in enumerate(screed.open(filename)):
                # make a new signature for each sequence
                sigs = signatures_factory()
                add_seq(sigs, record.sequence,
                        args.input_is_protein, args.check_sequence)

                set_sig_name(sigs, filename, name=record.name)
                siglist.extend(sigs)

            notify('calculated {} signatures for {} sequences in {}',
                   len(siglist), n + 1, filename)
        else:
            # make a single sig for the whole file
            sigs = signatures_factory()

            # consume & calculate signatures
            short_notify('... reading sequences from {}', filename)
            name = None
            for n, record in enumerate(screed.open(filename)):
                if n % 10000 == 0:
                    if n:
                        short_notify('\r...{} {}', filename, n, end='')
                    elif args.name_from_first:
                        name = record.name

                add_seq(sigs, record.sequence,
                        args.input_is_protein, args.check_sequence)

            short_notify('...{} {} sequences', filename, n, end='')

            set_sig_name(sigs, filename, name)
            siglist.extend(sigs)

            short_notify(f'calculated {len(siglist)} signatures for {n+1} sequences in {filename}')

        # if no --output specified, save to individual files w/in for loop
        if not args.output:
            save_siglist(siglist, sigfile)
            siglist = []

    # if --output specified, all collected signatures => args.output
    if args.output:
        save_siglist(siglist, args.output)
        siglist = []

    assert not siglist                    # juuuust checking.


def _compute_merged(args, signatures_factory):
    # make a signature for the whole file
    sigs = signatures_factory()

    n = 0
    total_seq = 0
    for filename in args.filenames:
        # consume & calculate signatures
        notify('... reading sequences from {}', filename)

        for n, record in enumerate(screed.open(filename)):
            if n % 10000 == 0 and n:
                notify('\r... {} {}', filename, n, end='')

            add_seq(sigs, record.sequence,
                    args.input_is_protein, args.check_sequence)
        notify('... {} {} sequences', filename, n + 1)

        total_seq += n + 1

    set_sig_name(sigs, filename, name=args.merge)
    notify('calculated 1 signature for {} sequences taken from {} files',
           total_seq, len(args.filenames))

    # at end, save!
    save_siglist(sigs, args.output)


def add_seq(sigs, seq, input_is_protein, check_sequence):
    for sig in sigs:
        if input_is_protein:
            sig.add_protein(seq)
        else:
            sig.add_sequence(seq, not check_sequence)


def set_sig_name(sigs, filename, name=None):
    for sig in sigs:
        if name is not None:
            sig._name = name
        sig.filename = filename


def save_siglist(siglist, sigfile_name):
    # save!
    with sourmash_args.FileOutput(sigfile_name, 'w') as fp:
        save_signatures(siglist, fp)
    notify('saved signature(s) to {}. Note: signature license is CC0.',
           sigfile_name)


class ComputeParameters(RustObject):
    __dealloc_func__ = lib.computeparams_free

    def __init__(self, ksizes, seed, protein, dayhoff, hp, dna, num_hashes, track_abundance, scaled):
        self._objptr = lib.computeparams_new()

        self.seed = seed
        self.ksizes = ksizes
        self.protein = protein
        self.dayhoff = dayhoff
        self.hp = hp
        self.dna = dna
        self.num_hashes = num_hashes
        self.track_abundance = track_abundance
        self.scaled = scaled

    @staticmethod
    def from_args(args):
        ptr = lib.computeparams_new()
        ret = ComputeParameters._from_objptr(ptr)

        for arg, value in vars(args).items():
            try:
                getattr(type(ret), arg).fset(ret, value)
            except AttributeError:
                pass

        return ret

    @property
    def seed(self):
        return self._methodcall(lib.computeparams_seed)

    @seed.setter
    def seed(self, v):
        return self._methodcall(lib.computeparams_set_seed, v)

    @property
    def ksizes(self):
        size = ffi.new("uintptr_t *")
        ksizes_ptr = self._methodcall(lib.computeparams_ksizes, size)
        size = size[0]
        ksizes = ffi.unpack(ksizes_ptr, size)
        lib.computeparams_ksizes_free(ksizes_ptr, size)
        return ksizes

    @ksizes.setter
    def ksizes(self, v):
        return self._methodcall(lib.computeparams_set_ksizes, list(v), len(v))

    @property
    def protein(self):
        return self._methodcall(lib.computeparams_protein)

    @protein.setter
    def protein(self, v):
        return self._methodcall(lib.computeparams_set_protein, v)

    @property
    def dayhoff(self):
        return self._methodcall(lib.computeparams_dayhoff)

    @dayhoff.setter
    def dayhoff(self, v):
        return self._methodcall(lib.computeparams_set_dayhoff, v)

    @property
    def hp(self):
        return self._methodcall(lib.computeparams_hp)

    @hp.setter
    def hp(self, v):
        return self._methodcall(lib.computeparams_set_hp, v)

    @property
    def dna(self):
        return self._methodcall(lib.computeparams_dna)

    @dna.setter
    def dna(self, v):
        return self._methodcall(lib.computeparams_set_dna, v)

    @property
    def num_hashes(self):
        return self._methodcall(lib.computeparams_num_hashes)

    @num_hashes.setter
    def num_hashes(self, v):
        return self._methodcall(lib.computeparams_set_num_hashes, v)

    @property
    def track_abundance(self):
        return self._methodcall(lib.computeparams_track_abundance)

    @track_abundance.setter
    def track_abundance(self, v):
        return self._methodcall(lib.computeparams_set_track_abundance, v)

    @property
    def scaled(self):
        return self._methodcall(lib.computeparams_scaled)

    @scaled.setter
    def scaled(self, v):
        return self._methodcall(lib.computeparams_set_scaled, int(v))
