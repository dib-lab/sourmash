"""
Tests for sourmash compute command-line functionality.
"""
import os
import gzip
import shutil
import screed
import glob
import json
import csv
import pytest

import sourmash_tst_utils as utils

import sourmash
from sourmash import MinHash
from sourmash.sbt import SBT, Node
from sourmash.sbtmh import SigLeaf, load_sbt_index
from sourmash.command_compute import ComputeParameters
from sourmash.cli.compute import subparser
from sourmash.cli import SourmashParser

from sourmash import signature
from sourmash import VERSION


@utils.in_tempdir
def test_do_sourmash_compute(c):
    testdata1 = utils.get_test_data("short.fa")

    c.run_sourmash("compute", "-k", "31", testdata1)

    sigfile = c.output("short.fa.sig")
    assert os.path.exists(sigfile)

    sig = next(signature.load_signatures(sigfile))
    assert str(sig).endswith('short.fa')


@utils.in_tempdir
def test_do_sourmash_compute_outdir(c):
    testdata1 = utils.get_test_data("short.fa")

    c.run_sourmash("compute", "-k", "31", testdata1, "--outdir", c.location)

    sigfile = c.output("short.fa.sig")
    assert os.path.exists(sigfile)

    sig = next(signature.load_signatures(sigfile))
    assert str(sig).endswith("short.fa")


@utils.in_tempdir
def test_do_sourmash_compute_output_valid_file(c):
    """ Trigger bug #123 """
    testdata1 = utils.get_test_data("short.fa")
    testdata2 = utils.get_test_data("short2.fa")
    testdata3 = utils.get_test_data("short3.fa")
    sigfile = c.output("short.fa.sig")

    c.run_sourmash(
        "compute", "-k", "31", "-o", sigfile, testdata1, testdata2, testdata3
    )

    assert os.path.exists(sigfile)
    assert not c.last_result.out  # stdout should be empty

    # is it valid json?
    with open(sigfile, "r") as f:
        data = json.load(f)

    filesigs = [sig["filename"] for sig in data]
    assert all(testdata in filesigs for testdata in (testdata1, testdata2, testdata3))


@utils.in_tempdir
def test_do_sourmash_compute_output_stdout_valid(c):
    testdata1 = utils.get_test_data("short.fa")
    testdata2 = utils.get_test_data("short2.fa")
    testdata3 = utils.get_test_data("short3.fa")

    c.run_sourmash("compute", "-k", "31", "-o", "-", testdata1, testdata2, testdata3)

    # is it valid json?
    data = json.loads(c.last_result.out)

    filesigs = [sig["filename"] for sig in data]
    assert all(testdata in filesigs for testdata in (testdata1, testdata2, testdata3))


@utils.in_tempdir
def test_do_sourmash_compute_output_and_name_valid_file(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    testdata3 = utils.get_test_data('short3.fa')
    sigfile = c.output('short.fa.sig')

    c.run_sourmash('compute', '-k', '31', '-o', sigfile, '--merge', '"name"', testdata1, testdata2, testdata3)

    assert os.path.exists(sigfile)
    assert 'calculated 1 signature for 4 sequences taken from 3 files' in c.last_result.err

    # is it valid json?
    with open(sigfile, 'r') as f:
        data = json.load(f)

    assert len(data) == 1

    sigfile_merged = c.output('short.all.fa.sig')
    c.run_sourmash('compute', '-k', '31', '-o', sigfile_merged, '--merge', '"name"', testdata1, testdata2, testdata3)

    with open(sigfile_merged, 'r') as f:
        data_merged = json.load(f)

    assert data[0]['signatures'][0]['mins'] == data_merged[0]['signatures'][0]['mins']


@utils.in_tempdir
def test_do_sourmash_compute_output_and_name_valid_file_outdir(c):
    testdata1 = utils.get_test_data('short.fa')
    testdata2 = utils.get_test_data('short2.fa')
    testdata3 = utils.get_test_data('short3.fa')
    sigfile = os.path.join(c.location, 'short.fa.sig')

    with pytest.raises(ValueError) as exc:
        c.run_sourmash('compute', '-k', '31', '-o', sigfile,
                       '--merge', '"name"',
                       testdata1, testdata2, testdata3,
                       '--outdir', c.location)

    errmsg = c.last_result.err
    assert "ERROR: --outdir doesn't make sense with -o/--output" in errmsg


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
        assert sig.name.endswith('shortName')


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
        assert sig.name == 'foo'

        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '31', '--name', 'foo',
                                            testdata1, '-o', 'foo2.sig'],
                                           in_directory=location)

        sigfile2 = os.path.join(location, 'foo2.sig')
        assert os.path.exists(sigfile2)

        sig2 = next(signature.load_signatures(sigfile))
        assert sig2.name == 'foo'
        assert sig.name == sig2.name


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
        assert sig.name == 'firstname'


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
        assert len(ksizes) == 2


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
            assert 7 in ksizes
            assert 10 in ksizes
            assert len(ksizes) == 4


def test_do_sourmash_compute_multik_with_dayhoff():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,30',
                                            '--dayhoff', '--no-dna',
                                            testdata1],
                                           in_directory=location)
        assert 'Computing only Dayhoff-encoded protein (and not nucleotide) ' \
               'signatures.' in err
        outfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(outfile)

        with open(outfile, 'rt') as fp:
            sigdata = fp.read()
            siglist = list(signature.load_signatures(sigdata))
            assert len(siglist) == 2
            ksizes = set([ x.minhash.ksize for x in siglist ])
            assert 7 in ksizes
            assert 10 in ksizes
            assert all(x.minhash.dayhoff for x in siglist)
            assert len(ksizes) == 2


def test_do_sourmash_compute_multik_with_dayhoff_and_dna():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,30',
                                            '--dayhoff',
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
            assert 7 in ksizes
            assert 10 in ksizes
            assert sum(x.minhash.moltype == 'DNA' for x in siglist) == 2
            assert sum(x.minhash.moltype == 'dayhoff' for x in siglist) == 2
            assert len(ksizes) == 4


def test_do_sourmash_compute_multik_with_hp():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,30',
                                            '--hp', '--no-dna',
                                            testdata1],
                                           in_directory=location)
        assert 'Computing only hp-encoded protein (and not nucleotide) ' \
               'signatures.' in err
        outfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(outfile)

        with open(outfile, 'rt') as fp:
            sigdata = fp.read()
            siglist = list(signature.load_signatures(sigdata))
            assert len(siglist) == 2
            ksizes = set([ x.minhash.ksize for x in siglist ])
            assert 7 in ksizes
            assert 10 in ksizes
            assert all(x.minhash.hp for x in siglist)
            assert len(ksizes) == 2


def test_do_sourmash_compute_multik_with_hp_and_dna():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,30',
                                            '--hp',
                                            testdata1],
                                           in_directory=location)
        outfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(outfile)

        with open(outfile, 'rt') as fp:
            sigdata = fp.read()
            siglist = list(signature.load_signatures(sigdata))
            assert len(siglist) == 4
            ksizes = set([ x.minhash.ksize for x in siglist ])
            assert 7 in ksizes
            assert 10 in ksizes
            assert 21 in ksizes
            assert 30 in ksizes
            assert len(ksizes) == 4


def test_do_sourmash_compute_multik_with_dayhoff_dna_protein():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,30',
                                            '--dayhoff', '--protein',
                                            testdata1],
                                           in_directory=location)
        outfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(outfile)

        with open(outfile, 'rt') as fp:
            sigdata = fp.read()
            siglist = list(signature.load_signatures(sigdata))
            assert len(siglist) == 6
            ksizes = set([ x.minhash.ksize for x in siglist ])
            assert 21 in ksizes
            assert 30 in ksizes
            assert 7 in ksizes
            assert 10 in ksizes
            assert sum(x.minhash.moltype == 'DNA' for x in siglist) == 2
            assert sum(x.minhash.moltype == 'dayhoff' for x in siglist) == 2
            assert sum(x.minhash.moltype == 'protein' for x in siglist) == 2
            assert len(ksizes) == 4


def test_do_sourmash_compute_multik_with_dayhoff_hp_dna_protein():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('short.fa')
        status, out, err = utils.runscript('sourmash',
                                           ['compute', '-k', '21,30',
                                            '--dayhoff', '--hp', '--protein',
                                            testdata1],
                                           in_directory=location)
        outfile = os.path.join(location, 'short.fa.sig')
        assert os.path.exists(outfile)

        with open(outfile, 'rt') as fp:
            sigdata = fp.read()
            siglist = list(signature.load_signatures(sigdata))
            assert len(siglist) == 8
            ksizes = set([ x.minhash.ksize for x in siglist ])
            assert 7 in ksizes
            assert 10 in ksizes
            assert 21 in ksizes
            assert 30 in ksizes
            assert sum(x.minhash.moltype == 'DNA' for x in siglist) == 2
            assert sum(x.minhash.moltype == 'dayhoff' for x in siglist) == 2
            assert sum(x.minhash.moltype == 'hp' for x in siglist) == 2
            # 2 = dayhoff, 2 = hp = 4 protein
            assert sum(x.minhash.moltype == 'protein' for x in siglist) == 2
            assert len(ksizes) == 4


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
        assert 7 in ksizes
        assert 10 in ksizes
        assert len(ksizes) == 2


def test_do_sourmash_compute_multik_protein_input_bad_ksize():
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
        assert status != 0
        assert 'protein ksizes must be divisible by 3' in err


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
        assert 7 in ksizes
        assert 10 in ksizes
        assert len(ksizes) == 2


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
            assert 7 in ksizes
            assert 10 in ksizes
            assert len(ksizes) == 2


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
            assert 7 in ksizes
            assert 10 in ksizes
            assert len(ksizes) == 2

            moltype = set([ x.minhash.moltype == 'protein'
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

        scaled_vals = [ x.minhash.scaled for x in siglist ]
        assert len(scaled_vals) == 2
        assert set(scaled_vals) == { 1 }


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

        max_hashes = [ x.minhash._max_hash for x in siglist ]
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

        max_hashes = [ x.minhash._max_hash for x in siglist ]
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
        sig1_aa, sig2_aa = sorted(x, key=lambda x: x.name)

        x = list(signature.load_signatures(sig2))
        sig1_trans, sig2_trans = sorted(x, key=lambda x: x.name)

        name1 = sig1_aa.name.split()[0]
        assert name1 == 'NP_414543.1'
        name2 = sig2_aa.name.split()[0]
        assert name2 == 'NP_414544.1'
        name3 = sig1_trans.name.split()[0]
        assert name3 == 'gi|556503834:2801-3733'
        name4 = sig2_trans.name.split()[0]
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
    sig1, sig2 = sorted(x, key=lambda x: x.name)
    print(sig1.name)
    print(sig2.name)

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
    sig1, sig2 = sorted(x, key=lambda x: x.name)

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
        sig1_aa, sig2_aa = sorted(x, key=lambda x: x.name)

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
        sig1_trans, sig2_trans = sorted(x, key=lambda x: x.name)

        knowngood = utils.get_test_data('benchmark.prot.sig')
        good_trans = list(signature.load_signatures(knowngood))[0]

        assert sig2_trans.similarity(good_trans) == 1.0


def test_compute_parameters():
    args_list = ["compute", "-k", "21,31", "--singleton", "--protein", "--no-dna", "input_file"]

    parser = SourmashParser(prog='sourmash')
    subp = parser.add_subparsers(title="instruction", dest="cmd", metavar="cmd")
    subparser(subp)

    args = parser.parse_args(args_list)

    params = ComputeParameters.from_args(args)

    assert params.ksizes == [21, 31]
    assert params.protein == True
    assert params.dna == False
    assert params.seed == 42
    assert params.dayhoff == False
    assert params.hp == False
    assert params.num_hashes == 500
    assert params.scaled == 0
    assert params.track_abundance == False
