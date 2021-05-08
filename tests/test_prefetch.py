"""
Tests for `sourmash prefetch` command-line and API functionality.
"""
import os
import csv
import pytest

import sourmash_tst_utils as utils
import sourmash


def test_prefetch_basic(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0

    assert "WARNING: no output(s) specified! Nothing will be saved from this prefetch!" in c.last_result.err
    assert "selecting specified query k=31" in c.last_result.err
    assert "loaded query: NC_009665.1 Shewanella baltica... (k=31, DNA)" in c.last_result.err
    assert "all sketches will be downsampled to scaled=1000" in c.last_result.err

    assert "total of 2 matching signatures." in c.last_result.err
    assert "of 5177 distinct query hashes, 5177 were found in matches above threshold." in c.last_result.err
    assert "a total of 0 query hashes remain unmatched." in c.last_result.err


def test_prefetch_query_abund(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch w/abund query
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0

    assert "WARNING: no output(s) specified! Nothing will be saved from this prefetch!" in c.last_result.err
    assert "selecting specified query k=31" in c.last_result.err
    assert "loaded query: NC_009665.1 Shewanella baltica... (k=31, DNA)" in c.last_result.err
    assert "all sketches will be downsampled to scaled=1000" in c.last_result.err

    assert "total of 2 matching signatures." in c.last_result.err
    assert "of 5177 distinct query hashes, 5177 were found in matches above threshold." in c.last_result.err
    assert "a total of 0 query hashes remain unmatched." in c.last_result.err


def test_prefetch_subj_abund(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch w/abund signature.
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0

    assert "WARNING: no output(s) specified! Nothing will be saved from this prefetch!" in c.last_result.err
    assert "selecting specified query k=31" in c.last_result.err
    assert "loaded query: NC_009665.1 Shewanella baltica... (k=31, DNA)" in c.last_result.err
    assert "all sketches will be downsampled to scaled=1000" in c.last_result.err

    assert "total of 2 matching signatures." in c.last_result.err
    assert "of 5177 distinct query hashes, 5177 were found in matches above threshold." in c.last_result.err
    assert "a total of 0 query hashes remain unmatched." in c.last_result.err


def test_prefetch_csv_out(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with CSV output
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    csvout = c.output('out.csv')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '-o', csvout, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(csvout)

    expected_intersect_bp = [2529000, 5177000]
    with open(csvout, 'rt', newline="") as fp:
        r = csv.DictReader(fp)
        for (row, expected) in zip(r, expected_intersect_bp):
            assert int(row['intersect_bp']) == expected


def test_prefetch_matches(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with --save-matches
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    matches_out = c.output('matches.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '--save-matches', matches_out, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(matches_out)

    sigs = sourmash.load_file_as_index(matches_out)

    expected_matches = [sig63, sig47]
    for (match, expected) in zip(sigs.signatures(), expected_matches):
        ss = sourmash.load_one_signature(expected, ksize=31)
        assert match == ss


def test_prefetch_matches_to_dir(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with --save-matches to a directory
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    ss63 = sourmash.load_one_signature(sig63)
    ss47 = sourmash.load_one_signature(sig47)

    matches_out = c.output('matches_dir/')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '--save-matches', matches_out, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(matches_out)
    assert os.path.isdir(matches_out)

    sigs = sourmash.load_file_as_signatures(matches_out)

    match_sigs = list(sigs)
    assert ss63 in match_sigs
    assert ss47 in match_sigs
    assert len(match_sigs) == 2


def test_prefetch_matches_to_sig_gz(runtmp, linear_gather):
    c = runtmp

    import gzip

    # test a basic prefetch, with --save-matches to a sig.gz file
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    ss63 = sourmash.load_one_signature(sig63)
    ss47 = sourmash.load_one_signature(sig47)

    matches_out = c.output('matches.sig.gz')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '--save-matches', matches_out, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(matches_out)
    assert os.path.isfile(matches_out)

    with gzip.open(matches_out, "rt") as fp:
        # can we read this as a gz file?
        fp.read()

    sigs = sourmash.load_file_as_signatures(matches_out)

    match_sigs = list(sigs)
    assert ss63 in match_sigs
    assert ss47 in match_sigs
    assert len(match_sigs) == 2


def test_prefetch_matches_to_zip(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with --save-matches to a zipfile
    import zipfile

    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')
    ss63 = sourmash.load_one_signature(sig63)
    ss47 = sourmash.load_one_signature(sig47)

    matches_out = c.output('matches.zip')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '--save-matches', matches_out, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(matches_out)
    assert os.path.isfile(matches_out)

    with zipfile.ZipFile(matches_out, "r") as fp:
        # can we read this as a .zip file?
        for zi in fp.infolist():
            pass

    sigs = sourmash.load_file_as_signatures(matches_out)

    match_sigs = list(sigs)
    assert ss63 in match_sigs
    assert ss47 in match_sigs
    assert len(match_sigs) == 2


def test_prefetch_matching_hashes(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with --save-matches
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    matches_out = c.output('matches.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63,
                   '--save-matching-hashes', matches_out, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(matches_out)

    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63, ksize=31)
    matches = set(ss47.minhash.hashes) & set(ss63.minhash.hashes)

    intersect = ss47.minhash.copy_and_clear()
    intersect.add_many(matches)

    ss = sourmash.load_one_signature(matches_out)
    assert ss.minhash == intersect


def test_prefetch_nomatch_hashes(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch, with --save-matches
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    nomatch_out = c.output('unmatched_hashes.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2,
                   '--save-unmatched-hashes', nomatch_out, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert os.path.exists(nomatch_out)

    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63, ksize=31)

    remain = ss47.minhash
    remain.remove_many(ss63.minhash.hashes)

    ss = sourmash.load_one_signature(nomatch_out)
    assert ss.minhash == remain


def test_prefetch_no_num_query(runtmp, linear_gather):
    c = runtmp

    # can't do prefetch with num signatures for query
    sig47 = utils.get_test_data('num/47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig47,
                       linear_gather)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status != 0


def test_prefetch_no_num_subj(runtmp, linear_gather):
    c = runtmp

    # can't do prefetch with num signatures for query; no matches!
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('num/63.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('prefetch', '-k', '31', sig47, sig63, linear_gather)

    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status != 0
    assert "ERROR in prefetch: no compatible signatures in any databases?!" in c.last_result.err


def test_prefetch_db_fromfile(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    from_file = c.output('from-list.txt')

    with open(from_file, 'wt') as fp:
        print(sig63, file=fp)
        print(sig2, file=fp)
        print(sig47, file=fp)

    c.run_sourmash('prefetch', '-k', '31', sig47, linear_gather,
                   '--db-from-file', from_file)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0

    assert "WARNING: no output(s) specified! Nothing will be saved from this prefetch!" in c.last_result.err
    assert "selecting specified query k=31" in c.last_result.err
    assert "loaded query: NC_009665.1 Shewanella baltica... (k=31, DNA)" in c.last_result.err
    assert "all sketches will be downsampled to scaled=1000" in c.last_result.err

    assert "total of 2 matching signatures." in c.last_result.err
    assert "of 5177 distinct query hashes, 5177 were found in matches above threshold." in c.last_result.err
    assert "a total of 0 query hashes remain unmatched." in c.last_result.err


def test_prefetch_no_db(runtmp, linear_gather):
    c = runtmp

    # test a basic prefetch with no databases/signatures
    sig47 = utils.get_test_data('47.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('prefetch', '-k', '31', sig47, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status != 0
    assert "ERROR: no databases or signatures to search!?" in c.last_result.err


def test_prefetch_downsample_scaled(runtmp, linear_gather):
    c = runtmp

    # test --scaled
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                   '--scaled', '1e5', linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "downsampling query from scaled=1000 to 10000" in c.last_result.err


def test_prefetch_empty(runtmp, linear_gather):
    c = runtmp

    # test --scaled
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    with pytest.raises(ValueError):
        c.run_sourmash('prefetch', '-k', '31', sig47, sig63, sig2, sig47,
                       '--scaled', '1e9', linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status != 0
    assert "no query hashes!? exiting." in c.last_result.err


def test_prefetch_basic_many_sigs(runtmp, linear_gather):
    c = runtmp

    # test what happens with many (and duplicate) signatures
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    manysigs = [sig63, sig2, sig47] * 5

    c.run_sourmash('prefetch', '-k', '31', sig47, *manysigs, linear_gather)
    print(c.last_result.status)
    print(c.last_result.out)
    print(c.last_result.err)

    assert c.last_result.status == 0
    assert "total of 10 matching signatures so far." in c.last_result.err
    assert "total of 10 matching signatures." in c.last_result.err
    assert "of 5177 distinct query hashes, 5177 were found in matches above threshold." in c.last_result.err
    assert "a total of 0 query hashes remain unmatched." in c.last_result.err
