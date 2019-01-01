"""
Tests for the 'sourmash signature' command line.
"""
from __future__ import print_function, unicode_literals
import csv

import pytest

from . import sourmash_tst_utils as utils
import sourmash

## command line tests


def test_run_sourmash_signature_cmd():
    status, out, err = utils.runscript('sourmash', ['signature'], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


def test_run_sourmash_sig_cmd():
    status, out, err = utils.runscript('sourmash', ['sig'], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


@utils.in_tempdir
def test_sig_merge_1(c):
    # merge of 47 & 63 should be union of mins
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    sig47and63 = utils.get_test_data('47+63.fa.sig')
    c.run_sourmash('sig', 'merge', sig47, sig63)

    # stdout should be new signature
    out = c.last_result.out

    test_merge_sig = sourmash.load_one_signature(sig47and63)
    actual_merge_sig = sourmash.load_one_signature(out)

    print(test_merge_sig.minhash)
    print(actual_merge_sig.minhash)
    print(out)

    assert actual_merge_sig.minhash == test_merge_sig.minhash


@utils.in_tempdir
def test_sig_merge_1_ksize_moltype(c):
    # check ksize, moltype args
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    sig47and63 = utils.get_test_data('47+63.fa.sig')
    c.run_sourmash('sig', 'merge', sig47, sig63, '--dna', '-k', '31')

    # stdout should be new signature
    out = c.last_result.out

    test_merge_sig = sourmash.load_one_signature(sig47and63)
    actual_merge_sig = sourmash.load_one_signature(out)

    print(test_merge_sig.minhash)
    print(actual_merge_sig.minhash)
    print(out)

    assert actual_merge_sig.minhash == test_merge_sig.minhash


@utils.in_tempdir
def test_sig_merge_2(c):
    # merge of 47 with nothing should be 47
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'merge', sig47)

    # stdout should be new signature
    out = c.last_result.out

    test_merge_sig = sourmash.load_one_signature(sig47)
    actual_merge_sig = sourmash.load_one_signature(out)

    print(out)

    assert actual_merge_sig.minhash == test_merge_sig.minhash


@utils.in_tempdir
def test_sig_merge_3_abund_ab_ok(c):
    # merge of 47 and 63 with abund should work
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    sig63abund = utils.get_test_data('track_abund/63.fa.sig')

    c.run_sourmash('sig', 'merge', sig47abund, sig63abund)
    actual_merge_sig = sourmash.load_one_signature(c.last_result.out)
    # @CTB: should check that this merge did what we think it should do!


@utils.in_tempdir
def test_sig_merge_3_abund_ab(c):
    # merge of 47 with abund, with 63 without, should fail; and vice versa
    sig47 = utils.get_test_data('47.fa.sig')
    sig63abund = utils.get_test_data('track_abund/63.fa.sig')

    with pytest.raises(ValueError) as e:
        c.run_sourmash('sig', 'merge', sig47, sig63abund)

    print(c.last_result)
    assert 'incompatible signatures: track_abundance is False in first sig, True in second' in c.last_result.err


@utils.in_tempdir
def test_sig_merge_3_abund_ba(c):
    # merge of 47 with abund, with 63 without, should fail; and vice versa
    sig47 = utils.get_test_data('47.fa.sig')
    sig63abund = utils.get_test_data('track_abund/63.fa.sig')

    with pytest.raises(ValueError) as e:
        c.run_sourmash('sig', 'merge', sig63abund, sig47)

    print(c.last_result)
    assert 'incompatible signatures: track_abundance is True in first sig, False in second' in c.last_result.err


@utils.in_tempdir
def test_sig_intersect_1(c):
    # intersect of 47 and 63 should be intersection of mins
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    sig47and63 = utils.get_test_data('47+63-intersect.fa.sig')
    c.run_sourmash('sig', 'intersect', sig47, sig63)

    # stdout should be new signature
    out = c.last_result.out

    test_intersect_sig = sourmash.load_one_signature(sig47and63)
    actual_intersect_sig = sourmash.load_one_signature(out)

    print(test_intersect_sig.minhash)
    print(actual_intersect_sig.minhash)
    print(out)

    assert actual_intersect_sig.minhash == test_intersect_sig.minhash


@utils.in_tempdir
def test_sig_intersect_2(c):
    # intersect of 47 with abund and 63 with abund should be same
    # as without abund, i.e. intersect 'flattens'
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')
    sig47and63 = utils.get_test_data('47+63-intersect.fa.sig')
    c.run_sourmash('sig', 'intersect', sig47, sig63)

    # stdout should be new signature
    out = c.last_result.out

    test_intersect_sig = sourmash.load_one_signature(sig47and63)
    actual_intersect_sig = sourmash.load_one_signature(out)

    print(test_intersect_sig.minhash)
    print(actual_intersect_sig.minhash)
    print(out)

    assert actual_intersect_sig.minhash == test_intersect_sig.minhash


@utils.in_tempdir
def test_sig_subtract_1(c):
    # subtract of 63 from 47
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    c.run_sourmash('sig', 'subtract', sig47, sig63)

    # stdout should be new signature
    out = c.last_result.out

    test1_sig = sourmash.load_one_signature(sig47)
    test2_sig = sourmash.load_one_signature(sig63)
    actual_subtract_sig = sourmash.load_one_signature(out)

    mins = set(test1_sig.minhash.get_mins())
    mins -= set(test2_sig.minhash.get_mins())

    assert set(actual_subtract_sig.minhash.get_mins()) == set(mins)


@utils.in_tempdir
def test_sig_subtract_2(c):
    # subtract of 63 from 47 should fail if 47 has abund
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('sig', 'subtract', sig47, sig63)


@utils.in_tempdir
def test_sig_subtract_3(c):
    # subtract of 63 from 47 should fail if 63 has abund
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('sig', 'subtract', sig47, sig63)


@utils.in_tempdir
def test_sig_intersect_2(c):
    # intersect of 47 and nothing should be self
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'intersect', sig47)

    # stdout should be new signature
    out = c.last_result.out

    test_intersect_sig = sourmash.load_one_signature(sig47)
    actual_intersect_sig = sourmash.load_one_signature(out)

    print(test_intersect_sig.minhash)
    print(actual_intersect_sig.minhash)
    print(out)

    assert actual_intersect_sig.minhash == test_intersect_sig.minhash


@utils.in_tempdir
def test_sig_rename_1(c):
    # set new name for 47
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'rename', sig47, 'fiz bar')

    # stdout should be new signature
    out = c.last_result.out

    test_rename_sig = sourmash.load_one_signature(sig47)
    actual_rename_sig = sourmash.load_one_signature(out)

    print(test_rename_sig.minhash)
    print(actual_rename_sig.minhash)

    assert actual_rename_sig.minhash == test_rename_sig.minhash
    assert test_rename_sig.name() != actual_rename_sig.name()
    assert actual_rename_sig.name() == 'fiz bar'


@utils.in_tempdir
def test_sig_extract_1(c):
    # extract 47 from 47... :)
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'extract', sig47)

    # stdout should be new signature
    out = c.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    assert actual_extract_sig == test_extract_sig


@utils.in_tempdir
def test_sig_extract_2(c):
    # extract matches to 47's md5sum from among several
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    c.run_sourmash('sig', 'extract', sig47, sig63, '--md5', '09a0869')

    # stdout should be new signature
    out = c.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    print(test_extract_sig.minhash)
    print(actual_extract_sig.minhash)

    assert actual_extract_sig == test_extract_sig


@utils.in_tempdir
def test_sig_extract_3(c):
    # extract nothing (no md5 match)
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'extract', sig47, '--md5', 'FOO')

    # stdout should be empty.
    out = c.last_result.out
    assert not out


@utils.in_tempdir
def test_sig_extract_4(c):
    # extract matches to 47's name from among several signatures
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    c.run_sourmash('sig', 'extract', sig47, sig63, '--name', 'NC_009665.1')

    # stdout should be new signature
    out = c.last_result.out

    test_extract_sig = sourmash.load_one_signature(sig47)
    actual_extract_sig = sourmash.load_one_signature(out)

    print(test_extract_sig.minhash)
    print(actual_extract_sig.minhash)

    assert actual_extract_sig == test_extract_sig


@utils.in_tempdir
def test_sig_extract_5(c):
    # extract nothing (no name match)
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'extract', sig47, '--name', 'FOO')

    # stdout should be empty.
    out = c.last_result.out
    assert not out


@utils.in_tempdir
def test_sig_extract_6(c):
    # extract matches to several names from among several signatures
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    c.run_sourmash('sig', 'extract', sig47, sig63, '--name', 'Shewanella')

    # stdout should be new signature
    out = c.last_result.out

    siglist = sourmash.load_signatures(out)
    siglist = list(siglist)

    assert len(siglist) == 2


@utils.in_tempdir
def test_sig_flatten_1(c):
    # extract matches to several names from among several signatures
    sig47abund = utils.get_test_data('track_abund/47.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'flatten', sig47abund, '--name', 'Shewanella')

    # stdout should be new signature
    out = c.last_result.out

    siglist = sourmash.load_signatures(out)
    siglist = list(siglist)

    assert len(siglist) == 1

    test_flattened = sourmash.load_one_signature(sig47)
    assert test_flattened.minhash == siglist[0].minhash


@utils.in_tempdir
def test_sig_downsample_1_scaled(c):
    # downsample a scaled signature
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'downsample', '--scaled', '10000', sig47)

    # stdout should be new signature
    out = c.last_result.out

    test_downsample_sig = sourmash.load_one_signature(sig47)
    actual_downsample_sig = sourmash.load_one_signature(out)

    test_mh = test_downsample_sig.minhash.downsample_scaled(10000)

    assert actual_downsample_sig.minhash == test_mh


@utils.in_tempdir
def test_sig_downsample_1_scaled_empty(c):
    # downsample a scaled signature
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('sig', 'downsample', sig47)


@utils.in_tempdir
def test_sig_downsample_2_num(c):
    # downsample a num signature
    sigs11 = utils.get_test_data('genome-s11.fa.gz.sig')
    c.run_sourmash('sig', 'downsample', '--num', '500',
                   '-k', '21', '--dna', sigs11)

    # stdout should be new signature
    out = c.last_result.out

    test_downsample_sig = sourmash.load_one_signature(sigs11, ksize=21,
                                                      select_moltype='DNA')
    actual_downsample_sig = sourmash.load_one_signature(out)
    test_mh = test_downsample_sig.minhash.downsample_n(500)

    assert actual_downsample_sig.minhash == test_mh


@utils.in_tempdir
def test_sig_downsample_2_num_to_scaled(c):
    # downsample a num signature and convert it into a scaled sig
    sigs11 = utils.get_test_data('genome-s11.fa.gz.sig')
    c.run_sourmash('sig', 'downsample', '--scaled', '10000',
                   '-k', '21', '--dna', sigs11)

    # stdout should be new signature
    out = c.last_result.out

    test_downsample_sig = sourmash.load_one_signature(sigs11, ksize=21,
                                                      select_moltype='DNA')
    actual_downsample_sig = sourmash.load_one_signature(out)

    test_mins = test_downsample_sig.minhash.get_mins()
    actual_mins = actual_downsample_sig.minhash.get_mins()

    # select those mins that are beneath the new max hash...
    max_hash = actual_downsample_sig.minhash.max_hash
    test_mins_down = { k for k in test_mins if k < max_hash }
    assert test_mins == actual_mins


@utils.in_tempdir
def test_sig_downsample_2_num_to_scaled_fail(c):
    # downsample a num signature and FAIL to convert it into a scaled sig
    # because new scaled is too low
    sigs11 = utils.get_test_data('genome-s11.fa.gz.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('sig', 'downsample', '--scaled', '100',
                       '-k', '21', '--dna', sigs11)


@utils.in_tempdir
def test_sig_downsample_2_num_empty(c):
    # downsample a num signature
    sigs11 = utils.get_test_data('genome-s11.fa.gz.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('sig', 'downsample', '-k', '21', '--dna', sigs11)


@utils.in_tempdir
def test_sig_info_1(c):
    # get basic info on a signature
    sig47 = utils.get_test_data('47.fa.sig')
    c.run_sourmash('sig', 'info', sig47)

    expected_output = """\
signature: NC_009665.1 Shewanella baltica OS185, complete genome
source file: 47.fa
md5: 09a08691ce52952152f0e866a59f6261
k=31 molecule=DNA num=0 scaled=1000 seed=42 track_abundance=0
size: 5177
signature license: CC0
""".splitlines()
    for line in expected_output:
        assert line.strip() in c.last_result.out


@utils.in_tempdir
def test_sig_info_2(c):
    # get info in CSV spreadsheet
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    c.run_sourmash('sig', 'info', sig47, sig63, '--csv', 'out.csv')


    expected_md5 = ['09a08691ce52952152f0e866a59f6261',
                    '38729c6374925585db28916b82a6f513']

    with open(c.output('out.csv'), 'rt') as fp:
        r = csv.DictReader(fp)

        n = 0

        for row, md5 in zip(r, expected_md5):
            assert row['md5'] == md5
            n += 1

        assert n == 2


@utils.in_tempdir
def test_import_export_1(c):
    # check to make sure we can import what we've exported!
    inp = utils.get_test_data('genome-s11.fa.gz.sig')
    outp = c.output('export.json')

    c.run_sourmash('sig', 'export', inp, '-o', outp, '-k', '21', '--dna')
    c.run_sourmash('sig', 'import', outp)

    original = sourmash.load_one_signature(inp, ksize=21, select_moltype='DNA')
    roundtrip = sourmash.load_one_signature(c.last_result.out)

    assert original.minhash == roundtrip.minhash


@utils.in_tempdir
def test_import_export_2(c):
    # check to make sure we can import a mash JSON dump file.
    # NOTE: msh.json_dump file calculated like so:
    #   mash sketch -s 500 -k 21 ./tests/test-data/genome-s11.fa.gz
    #   mash info -d ./tests/test-data/genome-s11.fa.gz.msh > tests/test-data/genome-s11.fa.gz.msh.json_dump
    #
    sig1 = utils.get_test_data('genome-s11.fa.gz.sig')
    msh_sig = utils.get_test_data('genome-s11.fa.gz.msh.json_dump')

    c.run_sourmash('sig', 'import', msh_sig)
    imported = sourmash.load_one_signature(c.last_result.out)
    compare = sourmash.load_one_signature(sig1, ksize=21, select_moltype='DNA')

    assert imported.minhash == compare.minhash
