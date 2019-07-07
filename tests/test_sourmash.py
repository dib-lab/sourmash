"""
Tests for the 'sourmash' command line.
"""
from __future__ import print_function, unicode_literals
import os
import gzip
import shutil
import time
import screed
import glob
import json
import csv
import pytest

from . import sourmash_tst_utils as utils
import sourmash
from sourmash import MinHash
from sourmash.sbt import SBT, Node
from sourmash.sbtmh import SigLeaf, load_sbt_index
try:
    import matplotlib
    matplotlib.use('Agg')
except ImportError:
    pass

from sourmash import signature
from sourmash import VERSION

def test_run_sourmash():
    status, out, err = utils.runscript('sourmash', [], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


def test_run_sourmash_badcmd():
    status, out, err = utils.runscript('sourmash', ['foobarbaz'], fail_ok=True)
    assert status != 0                    # bad arg!
    assert "Unrecognized command" in err

def test_sourmash_info():
    status, out, err = utils.runscript('sourmash', ['info'], fail_ok=False)

    # no output to stdout
    assert not out
    assert "sourmash version" in err
    assert "loaded from path" in err
    assert VERSION in err


def test_sourmash_info_verbose():
    status, out, err = utils.runscript('sourmash', ['info', '-v'])

    # no output to stdout
    assert not out
    assert "khmer version" in err
    assert "screed version" in err
    assert "loaded from path" in err


def test_do_sourmash_compute():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1],
                                           in_directory=location)

        sigfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(sigfile)

        sig = next(signature.load_signatures(sigfile))
        assert sig.name().endswith('short.fa')


def test_do_sourmash_compute_output_valid_file():
    """ Trigger bug #123 """
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        testdata3 = utils.get_test_data('short3.fa')
        sigfile = os.path.join(location, 'short.fa.sig')

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', '-o', sigfile,
                                            testdata1,
                                            testdata2, testdata3],
                                           in_directory=location)

        assert os.path.exists(sigfile)
        assert not out # stdout should be empty

        # is it valid json?
        with open(sigfile, 'r') as f:
            data = json.load(f)

        filesigs = [sig['filename'] for sig in data]
        assert all(testdata in filesigs
                   for testdata in (testdata1, testdata2, testdata3))


def test_do_sourmash_compute_output_stdout_valid():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        testdata3 = utils.get_test_data('short3.fa')

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', '-o', '-',
                                            testdata1,
                                            testdata2, testdata3],
                                           in_directory=location)

        # is it valid json?
        data = json.loads(out)

        filesigs = [sig['filename'] for sig in data]
        assert all(testdata in filesigs
                   for testdata in (testdata1, testdata2, testdata3))


def test_do_sourmash_compute_output_and_name_valid_file():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        testdata3 = utils.get_test_data('short3.fa')
        sigfile = os.path.join(location, 'short.fa.sig')

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', '-o', sigfile,
                                            '--merge', '"name"',
                                            testdata1,
                                            testdata2, testdata3],
                                           in_directory=location)

        assert os.path.exists(sigfile)
        assert 'calculated 1 signatures for 4 sequences taken from 3 files' in err

        # is it valid json?
        with open(sigfile, 'r') as f:
            data = json.load(f)

        assert len(data) == 1

        all_testdata = " ".join([testdata1, testdata2, testdata3])
        sigfile_merged = os.path.join(location, 'short.all.fa.sig')
        cmd = "cat {} | sourmash compute -k 31 -o {} -".format(
                all_testdata, sigfile_merged)
        status, out, err = utils.run_shell_cmd(cmd, in_directory=location)

        with open(sigfile_merged, 'r') as f:
            data_merged = json.load(f)

        assert data[0]['signatures'][0]['mins'] == data_merged[0]['signatures'][0]['mins']


def test_do_sourmash_compute_singleton():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', '--singleton',
                                            testdata1],
                                           in_directory=location)

        sigfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(sigfile)

        sig = next(signature.load_signatures(sigfile))
        assert sig.name().endswith('shortName')


def test_do_sourmash_compute_10x():
    bamnostic = pytest.importorskip('bamnostic')

    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('10x-example')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31',
                                            '--input-is-10x',
                                            testdata1],
                                           in_directory=location)

        sigfile = os.path.join(location, '10x-example.sig')
        assert os.path.exists(sigfile)

        with open(sigfile) as f:
            data = json.load(f)

        barcode_signatures = [sig['name'] for sig in data]

        with open(utils.get_test_data('10x-example/barcodes.tsv')) as f:
            true_barcodes = set(x.strip() for x in f.readlines())

        # Ensure that every cell barcode in barcodes.tsv has a signature
        assert all(bc in true_barcodes for bc in barcode_signatures)



def test_do_sourmash_compute_name():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', '--merge', 'foo',
                                            testdata1, '-o', 'foo.sig'],
                                           in_directory=location)

        sigfile = os.path.join(location, 'foo.sig')
        assert os.path.exists(sigfile)

        sig = next(signature.load_signatures(sigfile))
        assert sig.name() == 'foo'

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', '--name', 'foo',
                                            testdata1, '-o', 'foo2.sig'],
                                           in_directory=location)

        sigfile2 = os.path.join(location, 'foo2.sig')
        assert os.path.exists(sigfile)

        sig2 = next(signature.load_signatures(sigfile))
        assert sig2.name() == 'foo'
        assert sig.name() == sig2.name()


def test_do_sourmash_compute_name_fail_no_output():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', '--merge', 'foo',
                                            testdata1],
                                           in_directory=location,
                                           fail_ok=True)
        assert status == -1


def test_do_sourmash_compute_merge_fail_no_output():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', '--merge', 'foo',
                                            testdata1],
                                           in_directory=location,
                                           fail_ok=True)
        assert status == -1

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', '--name', 'foo',
                                            testdata1],
                                           in_directory=location,
                                           fail_ok=True)
        assert status == -1


def test_do_sourmash_compute_name_from_first():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short3.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', '--name-from-first',
                                            testdata1],
                                           in_directory=location)

        sigfile = os.path.join(location, 'short3.fa.sig')
        assert os.path.exists(sigfile)

        sig = next(signature.load_signatures(sigfile))
        assert sig.name() == 'firstname'


def test_do_sourmash_compute_multik():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            testdata1],
                                           in_directory=location)
        outfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(outfile)

        siglist = list(signature.load_signatures(outfile))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 21 in ksizes
        assert 31 in ksizes


def test_do_sourmash_compute_multik_with_protein():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,30',
                                            '--protein',
                                            testdata1],
                                           in_directory=location)
        outfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(outfile)

        with open(outfile, 'rt') as fp:
            sigdata = fp.read()
            siglist = list(signature.load_signatures(sigdata))
            assert len(siglist) == 4
            ksizes = set([ x.minhash.ksize for x in siglist ])
            assert 21 in ksizes
            assert 30 in ksizes


def test_do_sourmash_compute_multik_with_nothing():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--no-protein', '--no-dna',
                                            testdata1],
                                           in_directory=location,
                                           fail_ok=True)
        outfile = os.path.join(location, 'short.fa.sig')
        assert not os.path.exists(outfile)


def test_do_sourmash_compute_multik_protein_bad_ksize():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '20,32',
                                            '--protein', '--no-dna',
                                            testdata1],
                                           in_directory=location,
                                           fail_ok=True)
        outfile = os.path.join(location, 'short.fa.sig')
        assert not os.path.exists(outfile)
        assert 'protein ksizes must be divisible by 3' in err


@utils.in_tempdir
def test_do_sourmash_compute_multik_only_protein(c):
    # check sourmash compute with only protein, no nucl
    testdata1 = utils.get_test_data('short.fa')
    c.run_sourmash('compute', '-k', '21,30',
                   '--protein', '--no-dna', testdata1)
    outfile = os.path.join(c.location, 'short.fa.sig')
    assert os.path.exists(outfile)

    with open(outfile, 'rt') as fp:
        sigdata = fp.read()
        siglist = list(signature.load_signatures(sigdata))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 21 in ksizes
        assert 30 in ksizes


def test_do_sourmash_compute_multik_protein_input_non_div3_ksize():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short-protein.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '20,32',
                                            '--protein', '--no-dna',
                                            '--input-is-protein',
                                            testdata1],
                                           in_directory=location,
                                           fail_ok=True)
        outfile = os.path.join(location, 'short-protein.fa.sig')
        assert os.path.exists(outfile)


@utils.in_tempdir
def test_do_sourmash_compute_multik_only_protein_no_rna(c):
    # test --no-rna as well (otherwise identical to previous test)
    testdata1 = utils.get_test_data('short.fa')

    c.run_sourmash('compute', '-k', '21,30',
                   '--protein', '--no-rna', testdata1)
    outfile = os.path.join(c.location, 'short.fa.sig')
    assert os.path.exists(outfile)

    with open(outfile, 'rt') as fp:
        sigdata = fp.read()
        siglist = list(signature.load_signatures(sigdata))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 21 in ksizes
        assert 30 in ksizes


def test_do_sourmash_compute_protein_bad_sequences():
    """Proper error handling when Ns in dna sequence"""
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.bad.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,30',
                                            '--protein', '--no-dna',
                                            testdata1],
                                           in_directory=location)
        outfile = os.path.join(location, 'short.bad.fa.sig')
        assert os.path.exists(outfile)

        with open(outfile, 'rt') as fp:
            sigdata = fp.read()
            siglist = list(signature.load_signatures(sigdata))
            assert len(siglist) == 2
            ksizes = set([ x.minhash.ksize for x in siglist ])
            assert 21 in ksizes
            assert 30 in ksizes


def test_do_sourmash_compute_multik_input_is_protein():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('ecoli.faa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,30',
                                            '--input-is-protein',
                                            testdata1],
                                           in_directory=location)
        outfile = os.path.join(location, 'ecoli.faa.sig')
        assert os.path.exists(outfile)

        with open(outfile, 'rt') as fp:
            sigdata = fp.read()
            siglist = list(signature.load_signatures(sigdata))
            assert len(siglist) == 2
            ksizes = set([ x.minhash.ksize for x in siglist ])
            assert 21 in ksizes
            assert 30 in ksizes

            moltype = set([ x.minhash.is_molecule_type('protein')
                            for x in siglist ])
            assert len(moltype) == 1
            assert True in moltype


def test_do_sourmash_compute_multik_outfile():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        outfile = os.path.join(location, 'FOO.xxx')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            testdata1, '-o', outfile],
                                           in_directory=location)
        assert os.path.exists(outfile)

        siglist = list(signature.load_signatures(outfile))
        assert len(siglist) == 2
        ksizes = set([ x.minhash.ksize for x in siglist ])
        assert 21 in ksizes
        assert 31 in ksizes


def test_do_sourmash_compute_with_scaled_1():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        outfile = os.path.join(location, 'FOO.xxx')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--scaled', '1',
                                            testdata1, '-o', outfile],
                                            in_directory=location)
        assert os.path.exists(outfile)

        siglist = list(signature.load_signatures(outfile))
        assert len(siglist) == 2

        max_hashes = [ x.minhash.max_hash for x in siglist ]
        assert len(max_hashes) == 2
        assert set(max_hashes) == { sourmash.MAX_HASH }


def test_do_sourmash_compute_with_scaled_2():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        outfile = os.path.join(location, 'FOO.xxx')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--scaled', '2',
                                            testdata1, '-o', outfile],
                                            in_directory=location)
        assert os.path.exists(outfile)

        siglist = list(signature.load_signatures(outfile))
        assert len(siglist) == 2

        max_hashes = [ x.minhash.max_hash for x in siglist ]
        assert len(max_hashes) == 2
        assert set(max_hashes) == set([ int(2**64 /2.) ])


def test_do_sourmash_compute_with_scaled():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        outfile = os.path.join(location, 'FOO.xxx')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--scaled', '100',
                                            testdata1, '-o', outfile],
                                            in_directory=location)
        assert os.path.exists(outfile)

        siglist = list(signature.load_signatures(outfile))
        assert len(siglist) == 2

        max_hashes = [ x.minhash.max_hash for x in siglist ]
        assert len(max_hashes) == 2
        assert set(max_hashes) == set([ int(2**64 /100.) ])


def test_do_sourmash_compute_with_bad_scaled():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        outfile = os.path.join(location, 'FOO.xxx')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--scaled', '-1',
                                            testdata1, '-o', outfile],
                                            in_directory=location,
                                            fail_ok=True)

        assert status != 0
        assert '--scaled value must be >= 1' in err

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--scaled', '1000.5',
                                            testdata1, '-o', outfile],
                                            in_directory=location,
                                            fail_ok=True)

        assert status != 0
        assert '--scaled value must be integer value' in err

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--scaled', '1e9',
                                            testdata1, '-o', outfile],
                                            in_directory=location)

        assert status == 0
        assert 'WARNING: scaled value is nonsensical!?' in err


def test_do_sourmash_compute_with_seed():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        outfile = os.path.join(location, 'FOO.xxx')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,31',
                                            '--seed', '43',
                                            testdata1, '-o', outfile],
                                            in_directory=location)
        assert os.path.exists(outfile)

        siglist = list(signature.load_signatures(outfile))
        assert len(siglist) == 2

        seeds = [ x.minhash.seed for x in siglist ]
        assert len(seeds) == 2
        assert set(seeds) == set([ 43 ])


def test_do_sourmash_check_protein_comparisons():
    # this test checks 2 x 2 protein comparisons with E. coli genes.
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('ecoli.faa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21',
                                            '--input-is-protein',
                                            '--singleton',
                                            testdata1],
                                           in_directory=location)
        sig1 = os.path.join(location, 'ecoli.faa.sig')
        assert os.path.exists(sig1)

        testdata2 = utils.get_test_data('ecoli.genes.fna')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21',
                                            '--protein', '--no-dna',
                                            '--singleton',
                                            testdata2],
                                           in_directory=location)
        sig2 = os.path.join(location, 'ecoli.genes.fna.sig')
        assert os.path.exists(sig2)

        # I'm not sure why load_signatures is randomizing order, but ok.
        x = list(signature.load_signatures(sig1))
        sig1_aa, sig2_aa = sorted(x, key=lambda x: x.name())

        x = list(signature.load_signatures(sig2))
        sig1_trans, sig2_trans = sorted(x, key=lambda x: x.name())

        name1 = sig1_aa.name().split()[0]
        assert name1 == 'NP_414543.1'
        name2 = sig2_aa.name().split()[0]
        assert name2 == 'NP_414544.1'
        name3 = sig1_trans.name().split()[0]
        assert name3 == 'gi|556503834:2801-3733'
        name4 = sig2_trans.name().split()[0]
        assert name4 == 'gi|556503834:337-2799'

        print(name1, name3, round(sig1_aa.similarity(sig1_trans), 3))
        print(name2, name3, round(sig2_aa.similarity(sig1_trans), 3))
        print(name1, name4, round(sig1_aa.similarity(sig2_trans), 3))
        print(name2, name4, round(sig2_aa.similarity(sig2_trans), 3))

        assert round(sig1_aa.similarity(sig1_trans), 3) == 0.0
        assert round(sig2_aa.similarity(sig1_trans), 3) == 0.166
        assert round(sig1_aa.similarity(sig2_trans), 3) == 0.174
        assert round(sig2_aa.similarity(sig2_trans), 3) == 0.0


@utils.in_tempdir
def test_do_sourmash_check_knowngood_dna_comparisons(c):
    # this test checks against a known good signature calculated
    # by utils/compute-dna-mh-another-way.py
    testdata1 = utils.get_test_data('ecoli.genes.fna')
    c.run_sourmash('compute', '-k', '21',
                   '--singleton', '--dna',
                   testdata1)
    sig1 = c.output('ecoli.genes.fna.sig')
    assert os.path.exists(sig1)

    x = list(signature.load_signatures(sig1))
    sig1, sig2 = sorted(x, key=lambda x: x.name())

    knowngood = utils.get_test_data('benchmark.dna.sig')
    good = list(signature.load_signatures(knowngood))[0]

    assert sig2.similarity(good) == 1.0


@utils.in_tempdir
def test_do_sourmash_check_knowngood_dna_comparisons_use_rna(c):
    # check the --rna flag; otherwise identical to previous test.
    testdata1 = utils.get_test_data('ecoli.genes.fna')
    c.run_sourmash('compute', '-k', '21', '--singleton', '--rna',
                   testdata1)
    sig1 = c.output('ecoli.genes.fna.sig')
    assert os.path.exists(sig1)

    x = list(signature.load_signatures(sig1))
    sig1, sig2 = sorted(x, key=lambda x: x.name())

    knowngood = utils.get_test_data('benchmark.dna.sig')
    good = list(signature.load_signatures(knowngood))[0]

    assert sig2.similarity(good) == 1.0


def test_do_sourmash_check_knowngood_input_protein_comparisons():
    # this test checks against a known good signature calculated
    # by utils/compute-input-prot-another-way.py
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('ecoli.faa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21',
                                            '--input-is-protein',
                                            '--singleton',
                                            testdata1],
                                           in_directory=location)
        sig1 = os.path.join(location, 'ecoli.faa.sig')
        assert os.path.exists(sig1)

        x = list(signature.load_signatures(sig1))
        sig1_aa, sig2_aa = sorted(x, key=lambda x: x.name())

        knowngood = utils.get_test_data('benchmark.input_prot.sig')
        good_aa = list(signature.load_signatures(knowngood))[0]

        assert sig1_aa.similarity(good_aa) == 1.0


def test_do_sourmash_check_knowngood_protein_comparisons():
    # this test checks against a known good signature calculated
    # by utils/compute-prot-mh-another-way.py
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('ecoli.genes.fna')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21',
                                            '--singleton', '--protein',
                                            '--no-dna',
                                            testdata1],
                                           in_directory=location)
        sig1 = os.path.join(location, 'ecoli.genes.fna.sig')
        assert os.path.exists(sig1)

        x = list(signature.load_signatures(sig1))
        sig1_trans, sig2_trans = sorted(x, key=lambda x: x.name())

        knowngood = utils.get_test_data('benchmark.prot.sig')
        good_trans = list(signature.load_signatures(knowngood))[0]

        assert sig2_trans.similarity(good_trans) == 1.0


@utils.in_tempdir
def test_do_basic_compare(c):
    # try doing a basic compare
    import numpy
    testsigs = utils.get_test_data('genome-s1*.sig')
    testsigs = glob.glob(testsigs)

    c.run_sourmash('compare', '-o', 'cmp', '-k', '21', '--dna', *testsigs)

    cmp_outfile = c.output('cmp')
    assert os.path.exists(cmp_outfile)
    cmp_out = numpy.load(cmp_outfile.encode('utf-8'))

    sigs = []
    for fn in testsigs:
        sigs.append(sourmash.load_one_signature(fn, ksize=21,
                                                    select_moltype='dna'))

    cmp_calc = numpy.zeros([len(sigs), len(sigs)])
    for i, si in enumerate(sigs):
        for j, sj in enumerate(sigs):
            cmp_calc[i][j] = si.similarity(sj)

        sigs = []
        for fn in testsigs:
            sigs.append(sourmash.load_one_signature(fn, ksize=21,
                                                        select_moltype='dna'))
    assert (cmp_out == cmp_calc).all()


@utils.in_tempdir
def test_do_basic_compare_using_rna_arg(c):
    # try doing a basic compare using --rna instead of --dna
    import numpy
    testsigs = utils.get_test_data('genome-s1*.sig')
    testsigs = glob.glob(testsigs)

    c.run_sourmash('compare', '-o', 'cmp', '-k', '21', '--rna', *testsigs)

    cmp_outfile = c.output('cmp')
    assert os.path.exists(cmp_outfile)
    cmp_out = numpy.load(cmp_outfile.encode('utf-8'))

    sigs = []
    for fn in testsigs:
        sigs.append(sourmash.load_one_signature(fn, ksize=21,
                                                    select_moltype='dna'))

    cmp_calc = numpy.zeros([len(sigs), len(sigs)])
    for i, si in enumerate(sigs):
        for j, sj in enumerate(sigs):
            cmp_calc[i][j] = si.similarity(sj)

    assert (cmp_out == cmp_calc).all()


def test_do_compare_quiet():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31',
                                            testdata1, testdata2],
                                           in_directory=location)


        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '--csv', 'xxx',
                                            '-q'],
                                           in_directory=location)
        assert not out
        assert not err


def test_do_traverse_directory_compare():
    import numpy
    with utils.TempDirectory() as location:
        status, out, err = utils.runscript('sourmash',
                                           ['compare', '--traverse-directory',
                                            '-k 21', '--dna', utils.get_test_data('compare')],
                                           in_directory=location)
        print(out)
        assert 'genome-s10.fa.gz' in out
        assert 'genome-s11.fa.gz' in out

def test_do_compare_output_csv():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1, testdata2],
                                           in_directory=location)


        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '--csv', 'xxx'],
                                           in_directory=location)

        with open(os.path.join(location, 'xxx')) as fp:
            lines = fp.readlines()
            assert len(lines) == 3
            assert lines[1:] == ['1.0,0.93\n', '0.93,1.0\n']


def test_do_compare_downsample():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '--scaled', '200',
                                            '-k', '31', testdata1],
                                           in_directory=location)


        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '--scaled', '100',
                                            '-k', '31', testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '--csv', 'xxx'],
                                           in_directory=location)

        print(status, out, err)
        assert 'downsampling to scaled value of 200' in err
        with open(os.path.join(location, 'xxx')) as fp:
            lines = fp.readlines()
            assert len(lines) == 3
            assert lines[1].startswith('1.0,0.6666')
            assert lines[2].startswith('0.6666')


def test_do_compare_output_multiple_k():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21', testdata1],
                                           in_directory=location)
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '--csv', 'xxx'],
                                           in_directory=location,
                                           fail_ok=True)

        print(status, out, err)

        assert status == -1
        assert 'multiple k-mer sizes loaded; please specify one' in err
        assert '(saw k-mer sizes 21, 31)' in err


def test_do_compare_output_multiple_moltype():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21', '--dna', testdata1],
                                           in_directory=location)
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21', '--protein', testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '--csv', 'xxx'],
                                           in_directory=location,
                                           fail_ok=True)

        assert status == -1
        assert 'multiple molecule types loaded;' in err


def test_do_plot_comparison():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1, testdata2],
                                           in_directory=location)


        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '-o', 'cmp'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash', ['plot', 'cmp'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, "cmp.dendro.png"))
        assert os.path.exists(os.path.join(location, "cmp.matrix.png"))


def test_do_plot_comparison_2():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1, testdata2],
                                           in_directory=location)


        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '-o', 'cmp'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['plot', 'cmp', '--pdf'],
                                           in_directory=location)
        assert os.path.exists(os.path.join(location, "cmp.dendro.pdf"))
        assert os.path.exists(os.path.join(location, "cmp.matrix.pdf"))


def test_do_plot_comparison_3():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1, testdata2],
                                           in_directory=location)


        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig', '-o', 'cmp'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['plot', 'cmp', '--labels'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, "cmp.dendro.png"))
        assert os.path.exists(os.path.join(location, "cmp.matrix.png"))


def test_do_plot_comparison_5_force():
    import numpy
    with utils.TempDirectory() as location:
        D = numpy.zeros([2,2])
        D[0, 0] = 5
        with open(os.path.join(location, 'cmp'), 'wb') as fp:
            numpy.save(fp, D)

        with open(os.path.join(location, 'cmp.labels.txt'), 'wt') as fp:
            fp.write("a\nb\n")

        status, out, err = utils.runscript('sourmash',
                                           ['plot', 'cmp', '--labels', '-f'],
                                           in_directory=location)
        print(status, out, err)
        assert status == 0


def test_do_plot_comparison_4_fail_not_distance():
    import numpy
    with utils.TempDirectory() as location:
        D = numpy.zeros([2,2])
        D[0, 0] = 5
        with open(os.path.join(location, 'cmp'), 'wb') as fp:
            numpy.save(fp, D)

        with open(os.path.join(location, 'cmp.labels.txt'), 'wt') as fp:
            fp.write("a\nb\n")

        status, out, err = utils.runscript('sourmash',
                                           ['plot', 'cmp', '--labels'],
                                           in_directory=location, fail_ok=True)
        print(status, out, err)
        assert status != 0


def test_plot_subsample_1():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
        testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
        testdata4 = utils.get_test_data('genome-s10+s11.sig')
        inp_sigs = [testdata1, testdata2, testdata3, testdata4]

        status, out, err = utils.runscript('sourmash',
                                           ['compare'] + inp_sigs + \
                                           ['-o', 'cmp', '-k', '21', '--dna'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['plot', 'cmp',
                                            '--subsample', '3'],
                                           in_directory=location)

        print(out)

        expected = """\
0\ts10+s11
1\tgenome-s12.fa.gz
2\tgenome-s10.fa.gz"""
        assert expected in out


def test_plot_subsample_2():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
        testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
        testdata4 = utils.get_test_data('genome-s10+s11.sig')
        inp_sigs = [testdata1, testdata2, testdata3, testdata4]

        status, out, err = utils.runscript('sourmash',
                                           ['compare'] + inp_sigs + \
                                           ['-o', 'cmp', '-k', '21', '--dna'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['plot', 'cmp',
                                            '--subsample', '3',
                                            '--subsample-seed=2'],
                                           in_directory=location)

        print(out)
        expected = """\
0\tgenome-s12.fa.gz
1\ts10+s11
2\tgenome-s11.fa.gz"""
        assert expected in out


def test_search_query_sig_does_not_exist():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short2.fa.sig',
                                            'short.fa.sig'],
                                           in_directory=location, fail_ok=True)

        print(status, out, err)
        assert status == -1
        assert 'Cannot open file' in err
        assert len(err.splitlines()) < 5


def test_search_subject_sig_does_not_exist():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location, fail_ok=True)

        print(status, out, err)
        assert status == -1
        assert 'Cannot open file' in err


def test_search_second_subject_sig_does_not_exist():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'short.fa.sig', 'short2.fa.sig'],
                                           in_directory=location, fail_ok=True)

        print(status, out, err)
        assert status == -1
        assert 'Cannot open file' in err


def test_search():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1, testdata2],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches' in out
        assert '93.0%' in out


def test_search_ignore_abundance():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31',
                                            '--track-abundance',
                                            testdata1, testdata2],
                                           in_directory=location)



        # Make sure there's different percent matches when using or
        # not using abundance
        status1, out1, err1 = utils.runscript('sourmash',
                                           ['search',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)
        print(status1, out1, err1)
        assert '1 matches' in out1
        assert '81.5%' in out1

        status2, out2, err2 = utils.runscript('sourmash',
                                           ['search',
                                            '--ignore-abundance',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)
        print(status2, out2, err2)
        assert '1 matches' in out2
        assert '93.0%' in out2

        # Make sure results are different!
        assert out1 != out2


def test_search_csv():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1, testdata2],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'short2.fa.sig', '-o', 'xxx.csv'],
                                           in_directory=location)
        print(status, out, err)

        csv_file = os.path.join(location, 'xxx.csv')

        with open(csv_file) as fp:
            reader = csv.DictReader(fp)
            row = next(reader)
            assert float(row['similarity']) == 0.93
            assert row['name'].endswith('short2.fa')
            assert row['filename'].endswith('short2.fa.sig')
            assert row['md5'] == '914591cd1130aa915fe0c0c63db8f19d'


@utils.in_tempdir
def test_search_lca_db(c):
    # can we do a 'sourmash search' on an LCA database?
    query = utils.get_test_data('47.fa.sig')
    lca_db = utils.get_test_data('lca/47+63.lca.json')

    c.run_sourmash('search', query, lca_db)
    print(c)
    assert 'NC_009665.1 Shewanella baltica OS185, complete genome' in str(c)


@utils.in_tempdir
def test_gather_lca_db(c):
    # can we do a 'sourmash gather' on an LCA database?
    query = utils.get_test_data('47+63.fa.sig')
    lca_db = utils.get_test_data('lca/47+63.lca.json')

    c.run_sourmash('gather', query, lca_db)
    print(c)
    assert 'NC_009665.1 Shewanella baltica OS185' in str(c.last_result.out)


@utils.in_tempdir
def test_gather_csv_output_filename_bug(c):
    # check a bug where the database filename in the output CSV was incorrect
    query = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db_1 = utils.get_test_data('lca/delmont-1.lca.json')
    lca_db_2 = utils.get_test_data('lca/delmont-2.lca.json')

    c.run_sourmash('gather', query, lca_db_1, lca_db_2, '-o', 'out.csv')
    with open(c.output('out.csv'), 'rt') as fp:
        r = csv.DictReader(fp)
        row = next(r)
        assert row['filename'] == lca_db_1


def test_compare_deduce_molecule():
    # deduce DNA vs protein from query, if it is unique
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '30',
                                            '--no-dna', '--protein',
                                            testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert 'min similarity in matrix: 0.91' in out


def test_compare_choose_molecule_dna():
    # choose molecule type
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '30',
                                            '--dna', '--protein',
                                            testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compare', '--dna', 'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert 'min similarity in matrix: 0.938' in out


def test_compare_choose_molecule_protein():
    # choose molecule type
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '30',
                                            '--dna', '--protein',
                                            testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compare', '--protein', 'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert 'min similarity in matrix: 0.91' in out


def test_compare_no_choose_molecule_fail():
    # choose molecule type
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '30',
                                            '--dna', '--protein',
                                            testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location,
                                           fail_ok=True)

        assert 'multiple molecule types loaded; please specify' in err
        assert status != 0


def test_compare_deduce_ksize():
    # deduce ksize, if it is unique
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '29',
                                            testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compare', 'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert 'min similarity in matrix: 0.938' in out


def test_search_deduce_molecule():
    # deduce DNA vs protein from query, if it is unique
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '30',
                                            '--no-dna', '--protein',
                                            testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches' in out
        assert '(k=30, protein)' in err


def test_search_deduce_ksize():
    # deduce ksize from query, if it is unique
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '23',
                                            testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches' in out
        assert 'k=23' in err


def test_do_sourmash_index_multik_fail():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1],
                                           in_directory=location)
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '32', testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location, fail_ok=True)

        print(status, out, err)
        assert status == -1


def test_do_sourmash_index_multimol_fail():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1],
                                           in_directory=location)
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '--protein', '-k', '30',
                                            testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location, fail_ok=True)

        print(status, out, err)
        assert status == -1


def test_do_sourmash_index_multinum_fail():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', '-n', '500', testdata1],
                                           in_directory=location)
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', '-n', '1000', testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location, fail_ok=True)

        print(status, out, err)
        assert status == -1
        assert 'trying to build an SBT with incompatible signatures.' in err


def test_do_sourmash_index_multiscaled_fail():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '--scaled', '10', testdata1],
                                           in_directory=location)
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '--scaled', '1', testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location, fail_ok=True)

        print(status, out, err)
        assert status == -1
        assert 'trying to build an SBT with incompatible signatures.' in err


@utils.in_tempdir
def test_do_sourmash_index_multiscaled_rescale(c):
    # test sourmash index --scaled
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    c.run_sourmash('compute', '--scaled', '10', testdata1)
    c.run_sourmash('compute', '--scaled', '1', testdata2)

    c.run_sourmash('index', '-k', '31', 'zzz',
                   '--scaled', '10',
                   'short.fa.sig',
                   'short2.fa.sig')

    print(c)
    assert c.last_result.status == 0


@utils.in_tempdir
def test_do_sourmash_index_multiscaled_rescale_fail(c):
    # test sourmash index --scaled with invalid rescaling (10 -> 5)
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')

    c.run_sourmash('compute', '--scaled', '10', testdata1)
    c.run_sourmash('compute', '--scaled', '1', testdata2)
    # this should fail: cannot go from a scaled value of 10 to 5

    with pytest.raises(ValueError) as e:
        c.run_sourmash('index', '-k', '31', 'zzz',
                       '--scaled', '5',
                       'short.fa.sig',
                       'short2.fa.sig')

    print(e.value)
    assert c.last_result.status == -1
    assert 'new scaled 5 is lower than current sample scaled 10' in c.last_result.err


def test_do_sourmash_sbt_search_output():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', 'zzz', '-k', '31',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'zzz', '-o', 'foo'],
                                           in_directory=location)
        outfile = open(os.path.join(location, 'foo'))
        output = outfile.read()
        print(output)
        assert 'short.fa' in output
        assert 'short2.fa' in output


# check against a bug in sbt search triggered by incorrect max Jaccard
# calculation.
def test_do_sourmash_sbt_search_check_bug():
    with utils.TempDirectory() as location:
        # mins: 431
        testdata1 = utils.get_test_data('sbt-search-bug/nano.sig')

        # mins: 6264
        testdata2 = utils.get_test_data('sbt-search-bug/bacteroides.sig')

        status, out, err = utils.runscript('sourmash',
                                           ['index', 'zzz', '-k', '31',
                                            testdata1, testdata2],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', testdata1, 'zzz'],
                                           in_directory=location)
        assert '1 matches:' in out

        tree = load_sbt_index(os.path.join(location, 'zzz.sbt.json'))
        assert tree._nodes[0].metadata['min_n_below'] == 431


def test_do_sourmash_sbt_search_empty_sig():
    with utils.TempDirectory() as location:
        # mins: 431
        testdata1 = utils.get_test_data('sbt-search-bug/nano.sig')

        # mins: 0
        testdata2 = utils.get_test_data('sbt-search-bug/empty.sig')

        status, out, err = utils.runscript('sourmash',
                                           ['index', 'zzz', '-k', '31',
                                            testdata1, testdata2],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', testdata1, 'zzz'],
                                           in_directory=location)
        assert '1 matches:' in out

        tree = load_sbt_index(os.path.join(location, 'zzz.sbt.json'))
        assert tree._nodes[0].metadata['min_n_below'] == 1


def test_do_sourmash_sbt_move_and_search_output():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', 'zzz', '-k', '31',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        print(out)

        with open(os.path.join(location, 'zzz.sbt.json')) as fp:
            d = json.load(fp)
            assert d['storage']['args']['path'] == '.sbt.zzz'

        newpath = os.path.join(location, 'subdir')
        os.mkdir(newpath)

        # move both JSON file and subdirectory.
        shutil.move(os.path.join(location, 'zzz.sbt.json'), newpath)
        shutil.move(os.path.join(location, '.sbt.zzz'), newpath)

        status, out, err = utils.runscript('sourmash',
                                           ['search', '../short.fa.sig',
                                            'zzz', '-o', 'foo'],
                                           in_directory=newpath)
        outfile = open(os.path.join(newpath, 'foo'))
        output = outfile.read()
        print(output)
        assert 'short.fa' in output
        assert 'short2.fa' in output


def test_search_deduce_ksize_and_select_appropriate():
    # deduce ksize from query and select correct signature from DB
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '24',
                                            testdata1],
                                           in_directory=location)
        # The DB contains signatres for multiple ksizes
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '23,24',
                                            testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches' in out
        assert 'k=24' in err


def test_search_deduce_ksize_not_unique():
    # deduce ksize from query, fail because it is not unique
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '23,25',
                                            testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location,
                                           fail_ok=True)
        print(status, out, err)
        assert status == -1
        assert '2 signatures matching ksize' in err


def test_search_deduce_ksize_vs_user_specified():
    # user specified ksize is not available
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '23',
                                            testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['search', '-k', '24',
                                            'short.fa.sig', 'short2.fa.sig'],
                                           in_directory=location,
                                           fail_ok=True)
        print(status, out, err)
        assert status == -1
        assert '0 signatures matching ksize' in err


def test_search_containment():
    # search with --containment in signatures
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'short2.fa.sig', '--containment'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches' in out
        assert '95.8%' in out


def test_search_containment_sbt():
    # search with --containment in an SBT
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short2.fa.sig'],
                                           in_directory=location)


        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))
        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'zzz', '--containment'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches' in out
        assert '95.8%' in out


def test_search_gzip():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        data = open(os.path.join(location, 'short.fa.sig'), 'rb').read()
        with gzip.open(os.path.join(location, 'zzz.gz'), 'wb') as fp:
            fp.write(data)

        data = open(os.path.join(location, 'short2.fa.sig'), 'rb').read()
        with gzip.open(os.path.join(location, 'yyy.gz'), 'wb') as fp:
            fp.write(data)

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'zzz.gz',
                                            'yyy.gz'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches' in out
        assert '93.0%' in out


def test_search_2():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        testdata3 = utils.get_test_data('short3.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            testdata3],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'short2.fa.sig', 'short3.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '2 matches' in out
        assert '93.0%' in out
        assert '89.6%' in out


def test_search_3():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        testdata3 = utils.get_test_data('short3.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            testdata3],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['search', '-n', '1',
                                            'short.fa.sig',
                                            'short2.fa.sig', 'short3.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '2 matches; showing first 1' in out


def test_search_4():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        testdata3 = utils.get_test_data('short3.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            testdata3],
                                           in_directory=location)



        status, out, err = utils.runscript('sourmash',
                                           ['search', '-n', '0',
                                            'short.fa.sig',
                                            'short2.fa.sig', 'short3.fa.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '2 matches:' in out
        assert 'short2.fa' in out
        assert 'short3.fa' in out


def test_search_metagenome():
    with utils.TempDirectory() as location:
        testdata_glob = utils.get_test_data('gather/GCF*.sig')
        testdata_sigs = glob.glob(testdata_glob)

        query_sig = utils.get_test_data('gather/combined.sig')

        cmd = ['index', 'gcf_all', '-k', '21']
        cmd.extend(testdata_sigs)

        status, out, err = utils.runscript('sourmash', cmd,
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'gcf_all.sbt.json'))

        cmd = 'search {} gcf_all -k 21'.format(query_sig)
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        print(out)
        print(err)

        assert ' 33.2%       NC_003198.1 Salmonella enterica subsp. enterica serovar T...' in out
        assert '12 matches; showing first 3:' in out


def test_search_metagenome_traverse():
    with utils.TempDirectory() as location:
        testdata_dir = utils.get_test_data('gather')

        query_sig = utils.get_test_data('gather/combined.sig')

        cmd = 'search {} {} -k 21 --traverse-directory'
        cmd = cmd.format(query_sig, testdata_dir)
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        print(out)
        print(err)

        assert ' 33.2%       NC_003198.1 Salmonella enterica subsp. enterica serovar T...' in out
        assert '13 matches; showing first 3:' in out


# explanation: you cannot downsample a scaled SBT to match a scaled
# signature, so make sure that when you try such a search, it fails!
# (you *can* downsample a signature to match an SBT.)
def test_search_metagenome_downsample():
    with utils.TempDirectory() as location:
        testdata_glob = utils.get_test_data('gather/GCF*.sig')
        testdata_sigs = glob.glob(testdata_glob)

        query_sig = utils.get_test_data('gather/combined.sig')

        cmd = ['index', 'gcf_all', '-k', '21']
        cmd.extend(testdata_sigs)

        status, out, err = utils.runscript('sourmash', cmd,
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'gcf_all.sbt.json'))

        cmd = 'search {} gcf_all -k 21 --scaled 100000'.format(query_sig)
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location, fail_ok=True)
        assert status == -1


        assert "for tree 'gcf_all', scaled value is smaller than query." in err
        assert 'tree scaled: 10000; query scaled: 100000. Cannot do similarity search.' in err


def test_search_metagenome_downsample_containment():
    with utils.TempDirectory() as location:
        testdata_glob = utils.get_test_data('gather/GCF*.sig')
        testdata_sigs = glob.glob(testdata_glob)

        query_sig = utils.get_test_data('gather/combined.sig')

        cmd = ['index', 'gcf_all', '-k', '21']
        cmd.extend(testdata_sigs)

        status, out, err = utils.runscript('sourmash', cmd,
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'gcf_all.sbt.json'))

        cmd = 'search {} gcf_all -k 21 --scaled 100000 --containment'
        cmd = cmd.format(query_sig)
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        print(out)
        print(err)

        assert ' 32.9%       NC_003198.1 Salmonella enterica subsp. enterica serovar T...' in out
        assert '12 matches; showing first 3:' in out


@utils.in_tempdir
def test_search_metagenome_downsample_index(c):
    # does same search as search_metagenome_downsample_containment but
    # rescales during indexing
    #
    # for now, this test should fail; we need to clean up some internal
    # stuff before we can properly implement this!
    #
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_sig = utils.get_test_data('gather/combined.sig')

    # downscale during indexing, rather than during search.
    c.run_sourmash('index', 'gcf_all', '-k', '21', '--scaled', '100000',
                   *testdata_sigs)

    assert os.path.exists(c.output('gcf_all.sbt.json'))

    c.run_sourmash('search', query_sig, 'gcf_all', '-k', '21',
                       '--containment')
    print(c)

    assert ' 32.9%       NC_003198.1 Salmonella enterica subsp. enterica serovar T...' in str(c)
    assert ' 29.7%       NC_003197.2 Salmonella enterica subsp. enterica serovar T...' in str(c)
    assert '12 matches; showing first 3:' in str(c)


def test_mash_csv_to_sig():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa.msh.dump')
        testdata2 = utils.get_test_data('short.fa')

        status, out, err = utils.runscript('sourmash', ['import_csv',
                                                        testdata1,
                                                        '-o', 'xxx.sig'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31',
                                            '-n', '970', testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['search', '-k', '31',
                                            'short.fa.sig', 'xxx.sig'],
                                           in_directory=location)
        print(status, out, err)
        assert '1 matches:' in out
        assert '100.0%       short.fa' in out


def test_do_sourmash_index_bad_args():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig',
                                            '--dna', '--protein'],
                                           in_directory=location, fail_ok=True)

        print(out, err)
        assert "cannot specify both --dna/--rna and --protein!" in err
        assert status != 0


def test_do_sourmash_sbt_search():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'zzz'],
                                           in_directory=location)
        print(out)

        assert 'short.fa' in out
        assert 'short2.fa' in out


def test_do_sourmash_sbt_search_wrong_ksize():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            '-k', '31,51'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', '-k', '51',
                                            'short.fa.sig', 'zzz'],
                                           in_directory=location,
                                           fail_ok=True)

        assert status == -1
        assert 'this is different from' in err


def test_do_sourmash_sbt_search_multiple():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz2',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz2.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'zzz', 'zzz2'],
                                           in_directory=location)
        print(out)

        assert 'short.fa' in out
        assert 'short2.fa' in out


def test_do_sourmash_sbt_search_and_sigs():
    # search an SBT and a signature at same time.
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'zzz', 'short2.fa.sig'],
                                           in_directory=location)
        print(out)

        assert 'short.fa' in out
        assert 'short2.fa' in out


def test_do_sourmash_sbt_search_downsample():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            '--scaled=10'],
                                           in_directory=location)

        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1,
                                            '--scaled=5', '-o', 'query.sig'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'query.sig', 'zzz'],
                                           in_directory=location)
        print(out)

        assert 'short.fa' in out
        assert 'short2.fa' in out


def test_do_sourmash_sbt_search_downsample_2():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
        testdata2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')

        sbtname = 'foo'

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', sbtname,
                                            testdata2],
                                           in_directory=location)
        assert status == 0

        status, out, err = utils.runscript('sourmash',
                                           ['search', testdata1, sbtname,
                                            '--scaled=100000',
                                            '--threshold=0.01'],
                                           in_directory=location, fail_ok=True)
        assert status == -1
        assert 'Cannot do similarity search.' in err


def test_do_sourmash_index_single():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'zzz'],
                                           in_directory=location)
        print(out)

        assert 'short.fa' in out


def test_do_sourmash_sbt_search_selectprot():
    # index should fail when run on signatures with multiple types
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')

        args = ['compute', testdata1, testdata2,
                '--protein', '--dna', '-k', '30']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        args = ['index', '-k', '31', 'zzz', 'short.fa.sig', 'short2.fa.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location, fail_ok=True)

        print(out)
        print(err)
        assert status != 0


def test_do_sourmash_sbt_search_dnaprotquery():
    # sbt_search should fail if non-single query sig given
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')

        args = ['compute', testdata1, testdata2,
                '--protein', '--dna', '-k', '30']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        args = ['index', 'zzz', '--protein', 'short.fa.sig',
                'short2.fa.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        args = ['search', 'short.fa.sig', 'zzz']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location, fail_ok=True)
        assert status != 0
        print(out)
        print(err)
        assert 'need exactly one' in err


def test_do_sourmash_index_traverse():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            '--traverse-dir', '.'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))
        assert 'loaded 2 sigs; saving SBT under' in err

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'zzz'],
                                           in_directory=location)
        print(out)

        assert 'short.fa' in out
        assert 'short2.fa' in out


def test_do_sourmash_index_sparseness():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            '--traverse-dir', '.',
                                            '--sparseness', '1.0'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))
        assert 'loaded 2 sigs; saving SBT under' in err

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'short.fa.sig',
                                            'zzz'],
                                           in_directory=location)
        print(out)

        assert len(glob.glob(os.path.join(location, '.sbt.zzz', '*'))) == 2
        assert not glob.glob(os.path.join(location, '.sbt.zzz', '*internal*'))

        assert 'short.fa' in out
        assert 'short2.fa' in out


def test_do_sourmash_sbt_combine():
    with utils.TempDirectory() as location:
        files = [utils.get_test_data(f) for f in utils.SIG_FILES]

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz'] + files,
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['sbt_combine', 'joined',
                                            'zzz.sbt.json', 'zzz.sbt.json'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'joined.sbt.json'))

        filename = os.path.splitext(os.path.basename(utils.SIG_FILES[0]))[0]

        status, out, err = utils.runscript('sourmash',
                                           ['search', files[0], 'zzz'],
                                           in_directory=location)
        print(out)

        # we get notification of signature loading, too - so notify + result.
        assert out.count(filename) == 1

        status, out, err = utils.runscript('sourmash',
                                           ['search', files[0], 'joined'],
                                           in_directory=location)
        print(out)

        assert out.count(filename) == 1


def test_do_sourmash_index_append():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        testdata3 = utils.get_test_data('short3.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2, testdata3],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig', 'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        sbt_name = os.path.join(location, 'zzz',)
        sig_loc = os.path.join(location, 'short3.fa.sig')
        status, out, err = utils.runscript('sourmash',
                                           ['search', sig_loc, sbt_name])
        print(out)

        assert 'short.fa' in out
        assert 'short2.fa' in out
        assert 'short3.fa' not in out

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', '--append',
                                            'zzz',
                                            'short3.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        sbt_name = os.path.join(location, 'zzz',)
        sig_loc = os.path.join(location, 'short3.fa.sig')
        status, out, err = utils.runscript('sourmash',
                                           ['search',
                                            '--threshold', '0.95',
                                            sig_loc, sbt_name])
        print(out)

        assert 'short.fa' not in out
        assert 'short2.fa' in out
        assert 'short3.fa' in out


def test_do_sourmash_sbt_search_otherdir():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'xxx/zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'xxx', 'zzz.sbt.json'))

        sbt_name = os.path.join(location, 'xxx', 'zzz',)
        sig_loc = os.path.join(location, 'short.fa.sig')
        status, out, err = utils.runscript('sourmash',
                                           ['search', sig_loc, sbt_name])
        print(out)

        assert 'short.fa' in out
        assert 'short2.fa' in out


def test_do_sourmash_sbt_search_scaled_vs_num_1():
    # should not work: scaled query against num tree
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata2,
                                            '--scaled', '1000'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        sbt_name = os.path.join(location, 'zzz',)
        sig_loc = os.path.join(location, 'short2.fa.sig')
        status, out, err = utils.runscript('sourmash',
                                           ['search', sig_loc, sbt_name],
                                           fail_ok=True)

        assert status == -1
        assert 'tree and query are incompatible for search' in err


def test_do_sourmash_sbt_search_scaled_vs_num_2():
    # should not work: num query against scaled tree
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata2,
                                            '--scaled', '1000'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        sbt_name = os.path.join(location, 'zzz',)
        sig_loc = os.path.join(location, 'short.fa.sig')
        status, out, err = utils.runscript('sourmash',
                                           ['search', sig_loc, sbt_name],
                                           fail_ok=True)

        assert status == -1
        assert 'tree and query are incompatible for search' in err


def test_do_sourmash_sbt_search_scaled_vs_num_3():
    # should not work: scaled query against num signature
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata2,
                                            '--scaled', '1000'],
                                           in_directory=location)

        sig_loc = os.path.join(location, 'short.fa.sig')
        sig_loc2 = os.path.join(location, 'short2.fa.sig')
        status, out, err = utils.runscript('sourmash',
                                           ['search', sig_loc, sig_loc2],
                                           fail_ok=True)

        assert status == -1
        assert 'incompatible - cannot compare' in err


def test_do_sourmash_sbt_search_scaled_vs_num_4():
    # should not work: num query against scaled signature
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata2,
                                            '--scaled', '1000'],
                                           in_directory=location)

        sig_loc = os.path.join(location, 'short.fa.sig')
        sig_loc2 = os.path.join(location, 'short2.fa.sig')

        status, out, err = utils.runscript('sourmash',
                                           ['search', sig_loc2, sig_loc],
                                           fail_ok=True)
        assert status == -1
        assert 'incompatible - cannot compare' in err


def test_do_sourmash_check_search_vs_actual_similarity():
    with utils.TempDirectory() as location:
        files = [utils.get_test_data(f) for f in utils.SIG_FILES]

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz'] + files,
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        filename = os.path.splitext(os.path.basename(utils.SIG_FILES[0]))[0]

        status, out, err = utils.runscript('sourmash',
                                           ['search', files[0], 'zzz'],
                                           in_directory=location)
        assert status == 0


def test_do_sourmash_sbt_search_bestonly():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', '--best-only',
                                            'short.fa.sig', 'zzz'],
                                           in_directory=location)
        print(out)

        assert 'short.fa' in out


def test_sbt_search_order_dependence():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('genome-s10.fa.gz')
        testdata2 = utils.get_test_data('genome-s11.fa.gz')
        testdata3 = utils.get_test_data('genome-s12.fa.gz')
        testdata4 = utils.get_test_data('genome-s10+s11.fa.gz')

        cmd = 'compute --scaled 10000 -k 21,31 {} {} {} {}'
        cmd = cmd.format(testdata1, testdata2, testdata3, testdata4)

        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        cmd = 'index -k 21 134 genome-s10+s11.fa.gz.sig genome-s11.fa.gz.sig genome-s12.fa.gz.sig'
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        cmd = 'search -k 21 genome-s11.fa.gz.sig 134 --best-only -k 21 --dna'
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        print(out)
        print(err)
        assert '100.0%' in out


def test_sbt_search_order_dependence_2():
    # *should* return the same result as test_sbt_search_order_dependence,
    # but does not due to a bug.
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('genome-s10.fa.gz')
        testdata2 = utils.get_test_data('genome-s11.fa.gz')
        testdata3 = utils.get_test_data('genome-s12.fa.gz')
        testdata4 = utils.get_test_data('genome-s10+s11.fa.gz')

        cmd = 'compute --scaled 10000 -k 21,31 {} {} {} {}'
        cmd = cmd.format(testdata1, testdata2, testdata3, testdata4)

        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        cmd = 'index -k 21 314 genome-s11.fa.gz.sig genome-s10+s11.fa.gz.sig genome-s12.fa.gz.sig'
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        cmd = 'search -k 21 genome-s11.fa.gz.sig 314 --best-only --dna'
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        print(out)
        print(err)
        assert '100.0%' in out


def test_compare_with_abundance_1():
    with utils.TempDirectory() as location:
        # create two signatures
        E1 = MinHash(ksize=5, n=5, is_protein=False,
                        track_abundance=True)
        E2 = MinHash(ksize=5, n=5, is_protein=False,
                        track_abundance=True)

        E1.add_sequence('ATGGA')
        E2.add_sequence('ATGGA')

        s1 = signature.SourmashSignature(E1, filename='e1', name='e1')
        s2 = signature.SourmashSignature(E2, filename='e2', name='e2')

        signature.save_signatures([s1],
                                  open(os.path.join(location, 'e1.sig'), 'w'))
        signature.save_signatures([s2],
                                  open(os.path.join(location, 'e2.sig'), 'w'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'e1.sig', 'e2.sig',
                                            '-k' ,'5'],
                                           in_directory=location)
        assert '100.0%' in out


def test_compare_with_abundance_2():
    with utils.TempDirectory() as location:
        # create two signatures
        E1 = MinHash(ksize=5, n=5, is_protein=False,
                        track_abundance=True)
        E2 = MinHash(ksize=5, n=5, is_protein=False,
                        track_abundance=True)

        E1.add_sequence('ATGGA')

        E1.add_sequence('ATGGA')
        E2.add_sequence('ATGGA')

        s1 = signature.SourmashSignature(E1, filename='e1', name='e1')
        s2 = signature.SourmashSignature(E2, filename='e2', name='e2')

        signature.save_signatures([s1],
                                  open(os.path.join(location, 'e1.sig'), 'w'))
        signature.save_signatures([s2],
                                  open(os.path.join(location, 'e2.sig'), 'w'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'e1.sig', 'e2.sig',
                                            '-k' ,'5'],
                                           in_directory=location)
        assert '100.0%' in out


def test_compare_with_abundance_3():
    with utils.TempDirectory() as location:
        # create two signatures
        E1 = MinHash(ksize=5, n=5, is_protein=False,
                        track_abundance=True)
        E2 = MinHash(ksize=5, n=5, is_protein=False,
                        track_abundance=True)

        E1.add_sequence('ATGGA')
        E1.add_sequence('GGACA')

        E1.add_sequence('ATGGA')
        E2.add_sequence('ATGGA')

        s1 = signature.SourmashSignature(E1, filename='e1', name='e1')
        s2 = signature.SourmashSignature(E2, filename='e2', name='e2')

        signature.save_signatures([s1],
                                  open(os.path.join(location, 'e1.sig'), 'w'))
        signature.save_signatures([s2],
                                  open(os.path.join(location, 'e2.sig'), 'w'))

        status, out, err = utils.runscript('sourmash',
                                           ['search', 'e1.sig', 'e2.sig',
                                            '-k' ,'5'],
                                           in_directory=location)
        assert '70.5%' in out


def test_gather():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            '--scaled', '10'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata2,
                                            '--scaled', '10',
                                            '-o', 'query.fa.sig'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['gather',
                                            'query.fa.sig', 'zzz', '-o',
                                            'foo.csv', '--threshold-bp=1'],
                                           in_directory=location)

        print(out)
        print(err)

        assert '0.9 kbp      100.0%  100.0%' in out


def test_gather_csv():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            '--scaled', '10'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata2,
                                            '--scaled', '10',
                                            '-o', 'query.fa.sig'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['gather',
                                            'query.fa.sig', 'zzz', '-o',
                                            'foo.csv', '--threshold-bp=1'],
                                           in_directory=location)

        print(out)
        print(err)

        csv_file = os.path.join(location, 'foo.csv')

        with open(csv_file) as fp:
            reader = csv.DictReader(fp)
            row = next(reader)
            print(row)
            assert float(row['intersect_bp']) == 910
            assert float(row['f_orig_query']) == 1.0
            assert float(row['f_unique_to_query']) == 1.0
            assert float(row['f_match']) == 1.0
            assert row['name'].endswith('short2.fa')
            assert row['filename'] == 'zzz'
            assert row['md5'] == 'c9d5a795eeaaf58e286fb299133e1938'


def test_gather_multiple_sbts():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            '--scaled', '10'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata2,
                                            '--scaled', '10',
                                            '-o', 'query.fa.sig'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', 'zzz', '-k', '31',
                                            'short.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['index', 'zzz2', '-k', '31',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['gather',
                                            'query.fa.sig', 'zzz', 'zzz2',
                                            '-o', 'foo.csv',
                                            '--threshold-bp=1'],
                                           in_directory=location)

        print(out)
        print(err)

        assert '0.9 kbp      100.0%  100.0%' in out


def test_gather_sbt_and_sigs():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1, testdata2,
                                            '--scaled', '10'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata2,
                                            '--scaled', '10',
                                            '-o', 'query.fa.sig'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['gather',
                                            'query.fa.sig', 'zzz', 'short2.fa.sig',
                                            '-o', 'foo.csv',
                                            '--threshold-bp=1'],
                                           in_directory=location)

        print(out)
        print(err)

        assert '0.9 kbp      100.0%  100.0%' in out


def test_gather_file_output():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            '--scaled', '10'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata2,
                                            '--scaled', '10',
                                            '-o', 'query.fa.sig'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', '-k', '31', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['gather',
                                            'query.fa.sig', 'zzz',
                                            '--threshold-bp=500',
                                            '-o', 'foo.out'],
                                           in_directory=location)

        print(out)
        print(err)
        assert '0.9 kbp      100.0%  100.0%' in out
        with open(os.path.join(location, 'foo.out')) as f:
            output = f.read()
            print((output,))
            assert '910,1.0,1.0' in output


def test_gather_nomatch():
    with utils.TempDirectory() as location:
        testdata_query = utils.get_test_data('gather/GCF_000006945.2_ASM694v2_genomic.fna.gz.sig')
        testdata_match = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

        cmd = 'gather {} {}'.format(testdata_query, testdata_match)
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        print(out)
        print(err)

        assert 'found 0 matches total' in out
        assert 'the recovered matches hit 0.0% of the query' in out


def test_gather_metagenome():
    with utils.TempDirectory() as location:
        testdata_glob = utils.get_test_data('gather/GCF*.sig')
        testdata_sigs = glob.glob(testdata_glob)

        query_sig = utils.get_test_data('gather/combined.sig')

        cmd = ['index', 'gcf_all', '-k', '21']
        cmd.extend(testdata_sigs)

        status, out, err = utils.runscript('sourmash', cmd,
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'gcf_all.sbt.json'))

        cmd = 'gather {} gcf_all -k 21'.format(query_sig)
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        print(out)
        print(err)

        assert 'found 12 matches total' in out
        assert 'the recovered matches hit 100.0% of the query' in out
        assert all(('4.9 Mbp       33.2%  100.0%' in out,
                'NC_003198.1 Salmonella enterica subsp...' in out))
        assert all(('4.7 Mbp        0.5%    1.5%' in out,
                'NC_011294.1 Salmonella enterica subsp...' in out))

def test_multigather_metagenome():
    with utils.TempDirectory() as location:
        testdata_glob = utils.get_test_data('gather/GCF*.sig')
        testdata_sigs = glob.glob(testdata_glob)

        query_sig = utils.get_test_data('gather/combined.sig')

        cmd = ['index', 'gcf_all', '-k', '21']
        cmd.extend(testdata_sigs)

        status, out, err = utils.runscript('sourmash', cmd,
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'gcf_all.sbt.json'))

        cmd = 'multigather --query {} --db gcf_all -k 21'.format(query_sig)
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        print(out)
        print(err)

        assert 'found 12 matches total' in out
        assert 'the recovered matches hit 100.0% of the query' in out
        assert all(('4.9 Mbp       33.2%  100.0%' in out,
                'NC_003198.1 Salmonella enterica subsp...' in out))
        assert all(('4.7 Mbp        0.5%    1.5%' in out,
                'NC_011294.1 Salmonella enterica subsp...' in out))

def test_gather_metagenome_traverse():
    with utils.TempDirectory() as location:
        # set up a directory $location/gather that contains
        # everything in the 'tests/test-data/gather' directory
        # *except* the query sequence, which is 'combined.sig'.
        testdata_dir = utils.get_test_data('gather')
        copy_testdata = os.path.join(location, 'somesigs')
        shutil.copytree(testdata_dir, copy_testdata)
        os.unlink(os.path.join(copy_testdata, 'combined.sig'))

        query_sig = utils.get_test_data('gather/combined.sig')

        # now, feed in the new directory --
        cmd = 'gather {} {} -k 21 --traverse-directory'
        cmd = cmd.format(query_sig, copy_testdata)
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        print(cmd)
        print(out)
        print(err)

        assert 'found 12 matches total' in out
        assert 'the recovered matches hit 100.0% of the query' in out
        assert all(('4.9 Mbp       33.2%  100.0%' in out,
                'NC_003198.1 Salmonella enterica subsp...' in out))
        assert all(('4.7 Mbp        0.5%    1.5%' in out,
                'NC_011294.1 Salmonella enterica subsp...' in out))


def test_gather_metagenome_output_unassigned():
    with utils.TempDirectory() as location:
        testdata_glob = utils.get_test_data('gather/GCF_000195995*g')
        testdata_sigs = glob.glob(testdata_glob)[0]

        query_sig = utils.get_test_data('gather/combined.sig')

        cmd = 'gather {} {} -k 21'.format(query_sig, testdata_sigs)
        cmd += ' --output-unassigned=unassigned.sig'
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        print(out)
        print(err)

        assert 'found 1 matches total' in out
        assert 'the recovered matches hit 33.2% of the query' in out
        assert all(('4.9 Mbp       33.2%  100.0%' in out,
                'NC_003198.1 Salmonella enterica subsp...' in out))

        # now examine unassigned
        testdata2_glob = utils.get_test_data('gather/GCF_000009505.1*.sig')
        testdata2_sigs = glob.glob(testdata2_glob)[0]

        cmd = 'gather {} {} {} -k 21'.format('unassigned.sig', testdata_sigs, testdata2_sigs)
        status, out, err = utils.runscript('sourmash', cmd.split(' '),
                                           in_directory=location)

        print(out)
        print(err)
        assert all(('1.3 Mbp       13.6%   28.2%' in out,
                'NC_011294.1' in out))


@utils.in_tempdir
def test_gather_metagenome_output_unassigned_nomatches(c):
    # test --output-unassigned when there are no matches
    query_sig = utils.get_test_data('2.fa.sig')
    against_sig = utils.get_test_data('47.fa.sig')

    c.run_sourmash('gather', query_sig, against_sig,
                   '--output-unassigned', 'foo.sig')

    print(c.last_result.out)
    assert 'found 0 matches total;' in c.last_result.out

    x = sourmash.load_one_signature(query_sig, ksize=31)
    y = sourmash.load_one_signature(c.output('foo.sig'))

    assert x.minhash == y.minhash


def test_gather_metagenome_downsample():
    with utils.TempDirectory() as location:
        testdata_glob = utils.get_test_data('gather/GCF*.sig')
        testdata_sigs = glob.glob(testdata_glob)

        query_sig = utils.get_test_data('gather/combined.sig')

        cmd = ['index', 'gcf_all', '-k', '21']
        cmd.extend(testdata_sigs)

        status, out, err = utils.runscript('sourmash', cmd,
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'gcf_all.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['gather', query_sig, 'gcf_all',
                                            '-k', '21', '--scaled', '100000'],
                                           in_directory=location)

        print(out)
        print(err)

        assert 'found 11 matches total' in out
        assert 'the recovered matches hit 100.0% of the query' in out
        assert all(('5.2 Mbp       32.9%  100.0%' in out,
                'NC_003198.1' in out))
        assert all(('4.1 Mbp        0.6%    2.4%' in out,
                    '4.1 Mbp        4.4%   17.1%' in out))


def test_gather_query_downsample():
    with utils.TempDirectory() as location:
        testdata_glob = utils.get_test_data('gather/GCF*.sig')
        testdata_sigs = glob.glob(testdata_glob)

        query_sig = utils.get_test_data('GCF_000006945.2-s500.sig')

        status, out, err = utils.runscript('sourmash',
                                           ['gather', '-k', '31',
                                             query_sig] + testdata_sigs,
                                           in_directory=location)

        print(out)
        print(err)

        assert 'loaded 12 signatures' in err
        assert all(('4.9 Mbp      100.0%  100.0%' in out,
                'NC_003197.2' in out))


def test_gather_save_matches():
    with utils.TempDirectory() as location:
        testdata_glob = utils.get_test_data('gather/GCF*.sig')
        testdata_sigs = glob.glob(testdata_glob)

        query_sig = utils.get_test_data('gather/combined.sig')

        cmd = ['index', 'gcf_all', '-k', '21']
        cmd.extend(testdata_sigs)

        status, out, err = utils.runscript('sourmash', cmd,
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'gcf_all.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['gather', query_sig, 'gcf_all',
                                            '-k', '21',
                                            '--save-matches', 'save.sigs'],
                                           in_directory=location)

        print(out)
        print(err)

        assert 'found 12 matches total' in out
        assert 'the recovered matches hit 100.0% of the query' in out
        assert os.path.exists(os.path.join(location, 'save.sigs'))


def test_gather_error_no_cardinality_query():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1, testdata2],
                                           in_directory=location)

        testdata3 = utils.get_test_data('short3.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata3],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['gather',
                                            'short3.fa.sig', 'zzz'],
                                           in_directory=location,
                                           fail_ok=True)
        assert status == -1
        assert "query signature needs to be created with --scaled" in err


def test_gather_deduce_ksize():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            '--scaled', '10', '-k', '23'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata2,
                                            '--scaled', '10', '-k', '23',
                                            '-o', 'query.fa.sig'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['gather', 'query.fa.sig', 'zzz',
                                            '--threshold-bp=1'],
                                           in_directory=location)

        print(out)
        print(err)

        assert '0.9 kbp      100.0%  100.0%' in out


def test_gather_deduce_moltype():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        testdata2 = utils.get_test_data('short2.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata1, testdata2,
                                            '--scaled', '10', '-k', '30',
                                            '--no-dna', '--protein'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['compute', testdata2,
                                            '--scaled', '10', '-k', '30',
                                            '--no-dna', '--protein',
                                            '-o', 'query.fa.sig'],
                                           in_directory=location)

        status, out, err = utils.runscript('sourmash',
                                           ['index', 'zzz',
                                            'short.fa.sig',
                                            'short2.fa.sig'],
                                           in_directory=location)

        assert os.path.exists(os.path.join(location, 'zzz.sbt.json'))

        status, out, err = utils.runscript('sourmash',
                                           ['gather', 'query.fa.sig', 'zzz',
                                            '--threshold-bp=1'],
                                           in_directory=location)

        print(out)
        print(err)

        assert '1.9 kbp      100.0%  100.0%' in out


def test_gather_abund_1_1():
    # nullgraph/make-reads.py -S 1 -r 200 -C 2 tests/test-data/genome-s10.fa.gz > r1.fa
    # nullgraph/make-reads.py -S 1 -r 200 -C 20 tests/test-data/genome-s10.fa.gz > r2.fa
    # nullgraph/make-reads.py -S 1 -r 200 -C 2 tests/test-data/genome-s11.fa.gz > r3.fa
    # ./sourmash compute -k 21 --scaled 1000 --merge=1-1 -o reads-s10-s11.sig r[13].fa --track-abundance
    # ./sourmash compute -k 21 --scaled 1000 --merge=10-1 -o reads-s10x10-s11.sig r[23].fa --track-abundance

    with utils.TempDirectory() as location:
        query = utils.get_test_data('gather-abund/reads-s10-s11.sig')
        against_list = ['genome-s10', 'genome-s11', 'genome-s12']
        against_list = [ 'gather-abund/' + i + '.fa.gz.sig' \
                         for i in against_list ]
        against_list = [ utils.get_test_data(i) for i in against_list ]

        status, out, err = utils.runscript('sourmash',
                                           ['gather', query] + against_list,
                                           in_directory=location)

        print(out)
        print(err)

        assert '49.6%   78.5%       1.8    tests/test-data/genome-s10.fa.gz' in out
        assert '50.4%   80.0%       1.9    tests/test-data/genome-s11.fa.gz' in out
        assert 'genome-s12.fa.gz' not in out


def test_gather_abund_10_1():
    # nullgraph/make-reads.py -S 1 -r 200 -C 2 tests/test-data/genome-s10.fa.gz > r1.fa
    # nullgraph/make-reads.py -S 1 -r 200 -C 20 tests/test-data/genome-s10.fa.gz > r2.fa
    # nullgraph/make-reads.py -S 1 -r 200 -C 2 tests/test-data/genome-s11.fa.gz > r3.fa
    # ./sourmash compute -k 21 --scaled 1000 --merge=1-1 -o reads-s10-s11.sig r[13].fa --track-abundance
    # ./sourmash compute -k 21 --scaled 1000 --merge=10-1 -o reads-s10x10-s11.sig r[23].fa --track-abundance

    with utils.TempDirectory() as location:
        query = utils.get_test_data('gather-abund/reads-s10x10-s11.sig')
        against_list = ['genome-s10', 'genome-s11', 'genome-s12']
        against_list = [ 'gather-abund/' + i + '.fa.gz.sig' \
                         for i in against_list ]
        against_list = [ utils.get_test_data(i) for i in against_list ]

        status, out, err = utils.runscript('sourmash',
                                           ['gather', query, '-o', 'xxx.csv'] \
                                            + against_list,
                                           in_directory=location)

        print(out)
        print(err)
        assert '91.0%  100.0%      14.5    tests/test-data/genome-s10.fa.gz' in out
        assert '9.0%   80.0%       1.9    tests/test-data/genome-s11.fa.gz' in out
        assert 'genome-s12.fa.gz' not in out

        # check the calculations behind the above output by looking into
        # the CSV.
        with open(os.path.join(location, 'xxx.csv'), 'rt') as fp:
            r = csv.DictReader(fp)

            overlaps = []
            f_weighted_list = []
            average_abunds = []

            for row in r:
                overlap = float(row['intersect_bp'])
                f_weighted = float(row['f_unique_weighted'])
                average_abund = float(row['average_abund'])

                overlaps.append(overlap)
                f_weighted_list.append(f_weighted)
                average_abunds.append(average_abund)

            weighted_calc = []
            for (overlap, average_abund) in zip(overlaps, average_abunds):
                prod = overlap*average_abund
                weighted_calc.append(prod)

            total_weighted = sum(weighted_calc)
            for prod, f_weighted in zip(weighted_calc, f_weighted_list):
                assert prod / total_weighted == f_weighted, (prod, f_weighted)


def test_gather_abund_10_1_ignore_abundance():
    # nullgraph/make-reads.py -S 1 -r 200 -C 2 tests/test-data/genome-s10.fa.gz > r1.fa
    # nullgraph/make-reads.py -S 1 -r 200 -C 20 tests/test-data/genome-s10.fa.gz > r2.fa
    # nullgraph/make-reads.py -S 1 -r 200 -C 2 tests/test-data/genome-s11.fa.gz > r3.fa
    # ./sourmash compute -k 21 --scaled 1000 --merge=1-1 -o reads-s10-s11.sig r[13].fa --track-abundance
    # ./sourmash compute -k 21 --scaled 1000 --merge=10-1 -o reads-s10x10-s11.sig r[23].fa --track-abundance

    with utils.TempDirectory() as location:
        query = utils.get_test_data('gather-abund/reads-s10x10-s11.sig')
        against_list = ['genome-s10', 'genome-s11', 'genome-s12']
        against_list = [ 'gather-abund/' + i + '.fa.gz.sig' \
                         for i in against_list ]
        against_list = [ utils.get_test_data(i) for i in against_list ]

        status, out, err = utils.runscript('sourmash',
                                           ['gather', query] \
                                            + ['--ignore-abundance'] + \
                                            against_list,
                                           in_directory=location)

        print(out)
        print(err)
        assert all(('57.2%  100.0%', 'tests/test-data/genome-s10.fa.gz' in out))
        assert all(('42.8%   80.0%', 'tests/test-data/genome-s11.fa.gz' in out))
        assert 'genome-s12.fa.gz' not in out


def test_sbt_categorize():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
        testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
        testdata4 = utils.get_test_data('genome-s10+s11.sig')

        shutil.copyfile(testdata1, os.path.join(location, '1.sig'))
        shutil.copyfile(testdata2, os.path.join(location, '2.sig'))
        shutil.copyfile(testdata3, os.path.join(location, '3.sig'))
        shutil.copyfile(testdata4, os.path.join(location, '4.sig'))

        # omit 3
        args = ['index', '--dna', '-k', '21', 'zzz', '1.sig', '2.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        args = ['categorize', 'zzz', '--traverse-directory', '.',
                '--ksize', '21', '--dna', '--csv', 'out.csv']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        print(out)
        print(err)

        # mash dist genome-s10.fa.gz genome-s10+s11.fa.gz
        # yields 521/1000 ==> ~0.5
        assert 'for s10+s11, found: 0.50 genome-s10.fa.gz' in err

        out_csv = open(os.path.join(location, 'out.csv')).read()
        assert './4.sig,s10+s11,genome-s10.fa.gz,0.50' in out_csv


def test_sbt_categorize_ignore_abundance():
    with utils.TempDirectory() as location:

        query = utils.get_test_data('gather-abund/reads-s10x10-s11.sig')
        against_list = ['reads-s10-s11']
        against_list = [ 'gather-abund/' + i + '.sig' \
                         for i in against_list ]
        against_list = [ utils.get_test_data(i) for i in against_list ]

        # omit 3
        args = ['index', '--dna', '-k', '21', 'thebestdatabase'] + against_list
        status2, out2, err2 = utils.runscript('sourmash', args,
                                           in_directory=location)

        # --- Categorize without ignoring abundance ---
        args = ['categorize', 'thebestdatabase',
                '--ksize', '21', '--dna', '--csv', 'out3.csv', query]
        status3, out3, err3 = utils.runscript('sourmash', args,
                                           in_directory=location)

        print(out3)
        print(err3)

        assert 'for 1-1, found: 0.44 1-1' in err3

        out_csv3 = open(os.path.join(location, 'out3.csv')).read()
        assert 'reads-s10x10-s11.sig,1-1,1-1,0.4398' in out_csv3

        # --- Now categorize with ignored abundance ---
        args = ['categorize', '--ignore-abundance',
                '--ksize', '21', '--dna', '--csv', 'out4.csv',
                'thebestdatabase', query]
        status4, out4, err4 = utils.runscript('sourmash', args,
                                           in_directory=location)

        print(out4)
        print(err4)

        assert 'for 1-1, found: 0.88 1-1' in err4

        out_csv4 = open(os.path.join(location, 'out4.csv')).read()
        assert 'reads-s10x10-s11.sig,1-1,1-1,0.87699' in out_csv4

        # Make sure ignoring abundance produces a different output!
        assert err3 != err4


def test_sbt_categorize_already_done():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
        testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
        testdata4 = utils.get_test_data('genome-s10+s11.sig')

        shutil.copyfile(testdata1, os.path.join(location, '1.sig'))
        shutil.copyfile(testdata2, os.path.join(location, '2.sig'))
        shutil.copyfile(testdata3, os.path.join(location, '3.sig'))
        shutil.copyfile(testdata4, os.path.join(location, '4.sig'))

        # omit 3
        args = ['index', '--dna', '-k', '21', 'zzz', '1.sig', '2.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        with open(os.path.join(location, 'in.csv'), 'wt') as fp:
            fp.write('./4.sig,genome-s10.fa.gz,0.50')

        args = ['categorize', 'zzz', './2.sig', './4.sig',
                '--ksize', '21', '--dna', '--load-csv', 'in.csv']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        print(out)
        print(err)
        assert 'for genome-s11.fa.gz, no match found'
        assert not 'for s10+s11, found: 0.50 genome-s10.fa.gz' in err


def test_sbt_categorize_already_done_traverse():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
        testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')
        testdata4 = utils.get_test_data('genome-s10+s11.sig')

        shutil.copyfile(testdata1, os.path.join(location, '1.sig'))
        shutil.copyfile(testdata2, os.path.join(location, '2.sig'))
        shutil.copyfile(testdata3, os.path.join(location, '3.sig'))
        shutil.copyfile(testdata4, os.path.join(location, '4.sig'))

        # omit 3
        args = ['index', '--dna', '-k', '21', 'zzz', '1.sig', '2.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        with open(os.path.join(location, 'in.csv'), 'wt') as fp:
            fp.write('./4.sig,genome-s10.fa.gz,0.50')

        args = ['categorize', 'zzz', '--traverse-directory', '.',
                '--ksize', '21', '--dna', '--load-csv', 'in.csv']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        print(out)
        print(err)
        assert 'for genome-s11.fa.gz, no match found'
        assert not 'for s10+s11, found: 0.50 genome-s10.fa.gz' in err


def test_sbt_categorize_multiple_ksizes_moltypes():
    # 'categorize' should fail when there are multiple ksizes or moltypes
    # present
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        testdata2 = utils.get_test_data('genome-s11.fa.gz.sig')
        testdata3 = utils.get_test_data('genome-s12.fa.gz.sig')

        shutil.copyfile(testdata1, os.path.join(location, '1.sig'))
        shutil.copyfile(testdata2, os.path.join(location, '2.sig'))
        shutil.copyfile(testdata3, os.path.join(location, '3.sig'))

        args = ['index', '--dna', '-k', '21', 'zzz', '1.sig', '2.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        args = ['categorize', 'zzz', '--traverse-directory', '.']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location, fail_ok=True)

        assert status != 0
        assert 'multiple k-mer sizes/molecule types present' in err


def test_watch():
    with utils.TempDirectory() as location:
        testdata0 = utils.get_test_data('genome-s10.fa.gz')
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        shutil.copyfile(testdata1, os.path.join(location, '1.sig'))

        args = ['index', '--dna', '-k', '21', 'zzz', '1.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        cmd = """

             gunzip -c {} | sourmash watch --ksize 21 --dna zzz

        """.format(testdata0)
        status, out, err = utils.run_shell_cmd(cmd, in_directory=location)

        print(out)
        print(err)
        assert 'FOUND: genome-s10.fa.gz, at 1.000' in out


def test_watch_deduce_ksize():
    with utils.TempDirectory() as location:
        testdata0 = utils.get_test_data('genome-s10.fa.gz')
        utils.runscript('sourmash',
                        ['compute', testdata0, '-k', '29', '-o', '1.sig'],
                        in_directory=location)

        args = ['index', '--dna', '-k', '29', 'zzz', '1.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        cmd = """

             gunzip -c {} | sourmash watch --dna zzz

        """.format(testdata0)
        status, out, err = utils.run_shell_cmd(cmd, in_directory=location)

        print(out)
        print(err)
        assert 'Computing signature for k=29' in err
        assert 'genome-s10.fa.gz, at 1.000' in out


def test_watch_coverage():
    with utils.TempDirectory() as location:
        testdata0 = utils.get_test_data('genome-s10.fa.gz')
        testdata1 = utils.get_test_data('genome-s10.fa.gz.sig')
        shutil.copyfile(testdata1, os.path.join(location, '1.sig'))

        args = ['index', '--dna', '-k', '21', 'zzz', '1.sig']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        with open(os.path.join(location, 'query.fa'), 'wt') as fp:
            record = list(screed.open(testdata0))[0]
            for start in range(0, len(record), 100):
                fp.write('>{}\n{}\n'.format(start,
                                            record.sequence[start:start+500]))

        args = ['watch', '--ksize', '21', '--dna', 'zzz', 'query.fa']
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        print(out)
        print(err)
        assert 'FOUND: genome-s10.fa.gz, at 1.000' in out


def test_storage_convert_identity():
    with utils.TempDirectory() as location:
        testdata = utils.get_test_data('v2.sbt.json')
        shutil.copyfile(testdata, os.path.join(location, 'v2.sbt.json'))
        shutil.copytree(os.path.join(os.path.dirname(testdata), '.sbt.v2'),
                        os.path.join(location, '.sbt.v2'))
        testsbt = os.path.join(location, 'v2.sbt.json')

        original = SBT.load(testsbt, leaf_loader=SigLeaf.load)

        args = ['storage', 'convert', '-b', 'fsstorage', testsbt]
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        identity = SBT.load(testsbt, leaf_loader=SigLeaf.load)

        assert len(original) == len(identity)
        assert all(n1[1].name == n2[1].name
                   for (n1, n2) in zip(sorted(original), sorted(identity)))


def test_storage_convert_fsstorage_newpath():
    with utils.TempDirectory() as location:
        testdata = utils.get_test_data('v2.sbt.json')
        shutil.copyfile(testdata, os.path.join(location, 'v2.sbt.json'))
        shutil.copytree(os.path.join(os.path.dirname(testdata), '.sbt.v2'),
                        os.path.join(location, '.sbt.v2'))
        testsbt = os.path.join(location, 'v2.sbt.json')

        original = SBT.load(testsbt, leaf_loader=SigLeaf.load)

        args = ['storage', 'convert',
                           '-b', 'fsstorage({})'.format(os.path.join(location, 'v3')),
                           testsbt]
        status, out, err = utils.runscript('sourmash', args,
                                           in_directory=location)

        identity = SBT.load(testsbt, leaf_loader=SigLeaf.load)

        assert len(original) == len(identity)
        assert all(n1[1].name == n2[1].name
                   for (n1, n2) in zip(sorted(original), sorted(identity)))


def test_migrate():
    with utils.TempDirectory() as location:
        testdata = utils.get_test_data('v3.sbt.json')
        shutil.copyfile(testdata, os.path.join(location, 'v3.sbt.json'))
        shutil.copytree(os.path.join(os.path.dirname(testdata), '.sbt.v3'),
                        os.path.join(location, '.sbt.v3'))
        testsbt = os.path.join(location, 'v3.sbt.json')

        original = SBT.load(testsbt, leaf_loader=SigLeaf.load)

        status, out, err = utils.runscript('sourmash', ['migrate', testsbt],
                                           in_directory=location)

        identity = SBT.load(testsbt, leaf_loader=SigLeaf.load)

        assert len(original) == len(identity)
        assert all(n1[1].name == n2[1].name
                   for (n1, n2) in zip(sorted(original),
                                       sorted(identity)))

        assert "this is an old index version" not in err
        assert all('min_n_below' in node.metadata
                       for node in identity
                       if isinstance(node, Node))


def test_license_cc0():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', testdata1],
                                           in_directory=location)

        sigfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(sigfile)

        sig = next(signature.load_signatures(sigfile))
        assert sig.name().endswith('short.fa')

        assert sig.d['license'] == 'CC0'


def test_license_non_cc0():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31',
                                            testdata1, '--license', 'GPL'],
                                           in_directory=location, fail_ok=True)

        assert status != 0
        print(out)
        print(err)
        assert 'sourmash only supports CC0' in err


def test_license_load_non_cc0():
    with utils.TempDirectory() as location:
        sigfile = utils.get_test_data('bad-license.sig')

        try:
            sig = next(signature.load_signatures(sigfile, do_raise=True))
        except Exception as e:
            assert "sourmash only supports CC0-licensed signatures" in str(e)
