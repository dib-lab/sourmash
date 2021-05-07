"""
Tests for Index classes and subclasses.
"""
import pytest
import glob
import os
import zipfile
import shutil

import sourmash
from sourmash import load_one_signature, SourmashSignature
from sourmash.index import (LinearIndex, MultiIndex, ZipFileLinearIndex,
                            make_jaccard_search_query)
from sourmash.index.revindex import RevIndex
from sourmash.sbt import SBT, GraphFactory, Leaf
from sourmash.sbtmh import SigLeaf
from sourmash import sourmash_args
from sourmash.search import JaccardSearch, SearchType

import sourmash_tst_utils as utils


def test_simple_index(n_children):
    factory = GraphFactory(5, 100, 3)
    root = SBT(factory, d=n_children)

    leaf1_mh = sourmash.MinHash(0, 5, scaled=1)
    leaf1_mh.add_sequence("AAAAA")
    leaf1_mh.add_sequence("AAAAT")
    leaf1_mh.add_sequence("AAAAC")
    leaf1_sig = SourmashSignature(leaf1_mh)
    root.insert(leaf1_sig)

    leaf2_mh = sourmash.MinHash(0, 5, scaled=1)
    leaf2_mh.add_sequence("AAAAA")
    leaf2_mh.add_sequence("AAAAT")
    leaf2_mh.add_sequence("AAAAG")
    leaf2_sig = SourmashSignature(leaf2_mh)
    root.insert(leaf2_sig)
    
    leaf3_mh = sourmash.MinHash(0, 5, scaled=1)
    leaf3_mh.add_sequence("AAAAA")
    leaf3_mh.add_sequence("AAAAT")
    leaf3_mh.add_sequence("CAAAA")
    leaf3_sig = SourmashSignature(leaf3_mh)
    root.insert(leaf3_sig)
    
    leaf4_mh = sourmash.MinHash(0, 5, scaled=1)
    leaf4_mh.add_sequence("AAAAA")
    leaf4_mh.add_sequence("CAAAA")
    leaf4_mh.add_sequence("GAAAA")
    leaf4_sig = SourmashSignature(leaf4_mh)
    root.insert(leaf4_sig)
    
    leaf5_mh = sourmash.MinHash(0, 5, scaled=1)
    leaf5_mh.add_sequence("AAAAA")
    leaf5_mh.add_sequence("AAAAT")
    leaf5_mh.add_sequence("GAAAA")
    leaf5_sig = SourmashSignature(leaf5_mh)
    root.insert(leaf5_sig)
    
    linear = LinearIndex()
    linear.insert(leaf1_sig)
    linear.insert(leaf2_sig)
    linear.insert(leaf3_sig)
    linear.insert(leaf4_sig)
    linear.insert(leaf5_sig)

    search_fn = make_jaccard_search_query(do_containment=True)

    kmers = ["AAAAA", "AAAAT", "AAAAG", "CAAAA", "GAAAA"]
    for kmer in kmers:
        search_mh = sourmash.MinHash(0, 5, scaled=1)
        search_mh.add_sequence(kmer)
        search_sig = sourmash.SourmashSignature(search_mh)

        linear_found = linear.find(search_fn, search_sig)
        linear_found = set(linear_found)

        tree_found = set(root.find(search_fn, search_sig))
        
        assert tree_found
        assert tree_found == set(linear_found)


def test_linear_index_search():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = LinearIndex()
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    # now, search for sig2
    sr = lidx.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2

    # search for sig47 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss47, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss47
    assert sr[1][1] == ss63

    # search for sig63 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss63, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63
    assert sr[1][1] == ss47

    # search for sig63 with high threshold => 1 match
    sr = lidx.search(ss63, threshold=0.8)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63


def test_linear_index_gather():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = LinearIndex()
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    matches = lidx.gather(ss2)
    assert len(matches) == 1
    assert matches[0][0] == 1.0
    assert matches[0][1] == ss2

    matches = lidx.gather(ss47)
    assert len(matches) == 1
    assert matches[0][0] == 1.0
    assert matches[0][1] == ss47


def test_linear_index_search_subj_has_abundance():
    # check that signatures in the index are flattened appropriately.
    queryfile = utils.get_test_data('47.fa.sig')
    subjfile = utils.get_test_data('track_abund/47.fa.sig')

    qs = sourmash.load_one_signature(queryfile)
    ss = sourmash.load_one_signature(subjfile)

    linear = LinearIndex()
    linear.insert(ss)

    results = list(linear.search(qs, threshold=0))
    assert len(results) == 1
    # note: search returns _original_ signature, not flattened
    assert results[0].signature == ss


def test_linear_index_gather_subj_has_abundance():
    # check that signatures in the index are flattened appropriately.
    queryfile = utils.get_test_data('47.fa.sig')
    subjfile = utils.get_test_data('track_abund/47.fa.sig')

    qs = sourmash.load_one_signature(queryfile)
    ss = sourmash.load_one_signature(subjfile)

    linear = LinearIndex()
    linear.insert(ss)

    results = list(linear.gather(qs, threshold=0))
    assert len(results) == 1

    # note: gather returns _original_ signature, not flattened
    assert results[0].signature == ss


def test_index_search_subj_scaled_is_lower():
    # check that subject sketches are appropriately downsampled
    sigfile = utils.get_test_data('scaled100/GCF_000005845.2_ASM584v2_genomic.fna.gz.sig.gz')
    ss = sourmash.load_one_signature(sigfile)

    # double check :)
    assert ss.minhash.scaled == 100

    # build a new query that has a scaled of 1000
    qs = SourmashSignature(ss.minhash.downsample(scaled=1000))

    # create Index to search
    linear = LinearIndex()
    linear.insert(ss)

    # search!
    results = list(linear.search(qs, threshold=0))
    assert len(results) == 1
    # original signature (not downsampled) is returned
    assert results[0].signature == ss


def test_index_search_subj_num_is_lower():
    # check that subject sketches are appropriately downsampled
    sigfile = utils.get_test_data('num/47.fa.sig')
    ss = sourmash.load_one_signature(sigfile, ksize=31)

    # double check :)
    assert ss.minhash.num == 500

    # build a new query that has a num of 250
    qs = SourmashSignature(ss.minhash.downsample(num=250))

    # create Index to search
    linear = LinearIndex()
    linear.insert(ss)

    # search!
    results = list(linear.search(qs, threshold=0))
    assert len(results) == 1
    # original signature (not downsampled) is returned
    assert results[0].signature == ss


def test_index_search_query_num_is_lower():
    # check that query sketches are appropriately downsampled
    sigfile = utils.get_test_data('num/47.fa.sig')
    qs = sourmash.load_one_signature(sigfile, ksize=31)

    # double check :)
    assert qs.minhash.num == 500

    # build a new subject that has a num of 250
    ss = SourmashSignature(qs.minhash.downsample(num=250))

    # create Index to search
    linear = LinearIndex()
    linear.insert(ss)

    # search!
    results = list(linear.search(qs, threshold=0))
    assert len(results) == 1
    assert results[0].signature == ss


def test_linear_index_search_abund():
    # test Index.search_abund
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = LinearIndex()
    lidx.insert(ss47)
    lidx.insert(ss63)

    results = list(lidx.search_abund(ss47, threshold=0))
    assert len(results) == 2
    assert results[0].signature == ss47
    assert results[1].signature == ss63


def test_linear_index_search_abund_requires_threshold():
    # test Index.search_abund
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = LinearIndex()
    lidx.insert(ss47)
    lidx.insert(ss63)

    with pytest.raises(TypeError) as exc:
        results = list(lidx.search_abund(ss47, threshold=None))

    assert "'search_abund' requires 'threshold'" in str(exc.value)


def test_linear_index_search_abund_query_flat():
    # test Index.search_abund
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('track_abund/63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = LinearIndex()
    lidx.insert(ss47)
    lidx.insert(ss63)

    with pytest.raises(TypeError) as exc:
        results = list(lidx.search_abund(ss47, threshold=0))

    assert "'search_abund' requires query signature with abundance information" in str(exc.value)


def test_linear_index_search_abund_subj_flat():
    # test Index.search_abund
    sig47 = utils.get_test_data('track_abund/47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = LinearIndex()
    lidx.insert(ss47)
    lidx.insert(ss63)

    with pytest.raises(TypeError) as exc:
        results = list(lidx.search_abund(ss47, threshold=0))

    assert "'search_abund' requires subject signatures with abundance information" in str(exc.value)


def test_linear_index_save():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    linear = LinearIndex()
    linear.insert(ss2)
    linear.insert(ss47)
    linear.insert(ss63)
    
    with utils.TempDirectory() as location:
        filename = os.path.join(location, 'foo')
        linear.save(filename)

        si = set(sourmash.load_file_as_signatures(filename))

    x = {ss2, ss47, ss63}

    print(len(si))
    print(len(x))

    print(si)
    print(x)

    assert si == x, si


def test_linear_index_load():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    with utils.TempDirectory() as location:
        from sourmash import save_signatures

        filename = os.path.join(location, 'foo')
        with open(filename, 'wt') as fp:
            sourmash.save_signatures([ss2, ss47, ss63], fp)

        linear = LinearIndex.load(filename)

    x = {ss2, ss47, ss63}
    assert set(linear.signatures()) == x, linear.signatures
    assert linear.location == filename


def test_linear_index_save_load():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    linear = LinearIndex()
    linear.insert(ss2)
    linear.insert(ss47)
    linear.insert(ss63)
    
    with utils.TempDirectory() as location:
        filename = os.path.join(location, 'foo')
        linear.save(filename)
        linear2 = LinearIndex.load(filename)
        
    # now, search for sig2
    sr = linear2.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2


def test_linear_gather_threshold_1():
    # test gather() method, in some detail
    sig2 = load_one_signature(utils.get_test_data('2.fa.sig'), ksize=31)
    sig47 = load_one_signature(utils.get_test_data('47.fa.sig'), ksize=31)
    sig63 = load_one_signature(utils.get_test_data('63.fa.sig'), ksize=31)

    linear = LinearIndex()

    linear.insert(sig47)
    linear.insert(sig63)
    linear.insert(sig2)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(sig2.minhash.hashes.keys()))
    new_mh = sig2.minhash.copy_and_clear()

    # query with empty hashes
    assert not new_mh
    assert not linear.gather(SourmashSignature(new_mh))

    # add one hash
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 1

    results = linear.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig == sig2
    assert name is None

    # check with a threshold -> should be no results.
    results = linear.gather(SourmashSignature(new_mh), threshold_bp=5000)
    assert not results

    # add three more hashes => length of 4
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 4

    results = linear.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig == sig2
    assert name is None

    # check with a too-high threshold -> should be no results.
    results = linear.gather(SourmashSignature(new_mh), threshold_bp=5000)
    assert not results


def test_linear_gather_threshold_5():
    # test gather() method above threshold
    sig2 = load_one_signature(utils.get_test_data('2.fa.sig'), ksize=31)
    sig47 = load_one_signature(utils.get_test_data('47.fa.sig'), ksize=31)
    sig63 = load_one_signature(utils.get_test_data('63.fa.sig'), ksize=31)

    linear = LinearIndex(filename='foo')

    linear.insert(sig47)
    linear.insert(sig63)
    linear.insert(sig2)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(sig2.minhash.hashes.keys()))
    new_mh = sig2.minhash.copy_and_clear()

    # add five hashes
    for i in range(5):
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())
        new_mh.add_hash(mins.pop())

    # should get a result with no threshold (any match at all is returned)
    results = linear.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig == sig2
    assert name == 'foo'

    # now, check with a threshold_bp that should be meet-able.
    results = linear.gather(SourmashSignature(new_mh), threshold_bp=5000)
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig == sig2
    assert name == 'foo'


def test_linear_index_multik_select():
    # this loads three ksizes, 21/31/51
    sig2 = utils.get_test_data('2.fa.sig')
    siglist = sourmash.load_file_as_signatures(sig2)

    linear = LinearIndex()
    for ss in siglist:
        linear.insert(ss)

    # select most specifically
    linear2 = linear.select(ksize=31, moltype='DNA')
    assert len(linear2) == 1

    # all are DNA:
    linear2 = linear.select(moltype='DNA')
    assert len(linear2) == 3


def test_linear_index_moltype_select():
    # this loads two ksizes(21, 10), and two moltypes (DNA and protein)
    filename = utils.get_test_data('genome-s10+s11.sig')
    siglist = sourmash.load_file_as_signatures(filename)

    linear = LinearIndex()
    for ss in siglist:
        linear.insert(ss)

    # select most specific DNA
    linear2 = linear.select(ksize=30, moltype='DNA')
    assert len(linear2) == 1

    # select most specific protein
    linear2 = linear.select(ksize=10, moltype='protein')
    assert len(linear2) == 1

    # can leave off ksize, selects all ksizes
    linear2 = linear.select(moltype='DNA')
    assert len(linear2) == 2

    # can leave off ksize, selects all ksizes
    linear2 = linear.select(moltype='protein')
    assert len(linear2) == 2

    # select something impossible
    linear2 = linear.select(ksize=4)
    assert len(linear2) == 0


@utils.in_tempdir
def test_index_same_md5sum_fsstorage(c):
    testdata1 = utils.get_test_data('img/2706795855.sig')
    testdata2 = utils.get_test_data('img/638277004.sig')

    c.run_sourmash('index', '-k', '21', 'zzz.sbt.json', testdata1, testdata2)
    assert c.last_result.status == 0

    outfile = c.output('zzz.sbt.json')
    assert os.path.exists(outfile)
    storage = c.output('.sbt.zzz')
    assert len(glob.glob(storage + "/*")) == 3


@utils.in_tempdir
def test_index_same_md5sum_sbt_zipstorage(c):
    testdata1 = utils.get_test_data('img/2706795855.sig')
    testdata2 = utils.get_test_data('img/638277004.sig')

    c.run_sourmash('index', '-k', '21', 'zzz.sbt.zip', testdata1, testdata2)
    assert c.last_result.status == 0

    outfile = c.output('zzz.sbt.zip')
    assert os.path.exists(outfile)
    zout = zipfile.ZipFile(outfile, mode='r')
    # should have 3 files, 1 internal and two sigs. We check for 4 because the
    # directory also shows in namelist()
    assert len([f for f in zout.namelist() if f.startswith(".sbt.zzz/")]) == 4


@utils.in_thisdir
def test_zipfile_protein_command_search(c):
    # test command-line search/gather of zipfile with protein sigs
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/protein.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches:' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out)
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


@utils.in_thisdir
def test_zipfile_hp_command_search(c):
    # test command-line search/gather of zipfile with hp sigs
    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/hp.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches:' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


@utils.in_thisdir
def test_zipfile_dayhoff_command_search(c):
    # test command-line search/gather of zipfile with dayhoff sigs
    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/dayhoff.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches:' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


@utils.in_thisdir
def test_zipfile_protein_command_search_combined(c):
    # test command-line search/gather of combined zipfile with protein sigs
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/all.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches:' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out)
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


@utils.in_thisdir
def test_zipfile_hp_command_search_combined(c):
    # test command-line search/gather of combined zipfile with hp sigs
    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/all.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches:' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


@utils.in_thisdir
def test_zipfile_dayhoff_command_search_combined(c):
    # test command-line search/gather of combined zipfile with dayhoff sigs
    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/all.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches:' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


@utils.in_thisdir
def test_zipfile_dayhoff_command_search_protein(c):
    # test command-line search/gather of protein sigs in zipfile
    # with dayhoff query
    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/protein.zip')

    with pytest.raises(ValueError) as exc:
        c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')

    assert 'no compatible signatures found in ' in c.last_result.err


def test_zipfile_API_signatures():
    # return all of the .sig and .sig.gz files in all.zip
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db)
    siglist = list(zipidx.signatures())
    assert len(siglist) == 7
    assert len(zipidx) == 7


def test_zipfile_API_signatures_traverse_yield_all():
    # include dna-sig.noext, but not build.sh (cannot be loaded as signature)
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db, traverse_yield_all=True)
    siglist = list(zipidx.signatures())
    assert len(siglist) == 8
    assert len(zipidx) == 8

    # confirm that there are 12 files in there total, incl build.sh and dirs
    zf = zipidx.zf
    allfiles = [ zi.filename for zi in zf.infolist() ]
    print(allfiles)
    assert len(allfiles) == 12


def test_zipfile_API_signatures_traverse_yield_all_select():
    # include dna-sig.noext
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db, traverse_yield_all=True)
    zipidx = zipidx.select(moltype='DNA')
    siglist = list(zipidx.signatures())
    assert len(siglist) == 2
    assert len(zipidx) == 2


def test_zipfile_API_signatures_select():
    # include dna-sig.noext
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db)
    zipidx = zipidx.select(moltype='DNA')
    siglist = list(zipidx.signatures())
    assert len(siglist) == 1
    assert len(zipidx) == 1


def test_zipfile_API_save():
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db)

    with pytest.raises(NotImplementedError):
        zipidx.save('xxx')


def test_zipfile_API_insert():
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db)

    with pytest.raises(NotImplementedError):
        # at some point probably want to change this to a real signature :)
        zipidx.insert(None)


def test_zipfile_API_location():
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db)

    assert zipidx.location == zipfile_db


def test_zipfile_load_file_as_signatures():
    from types import GeneratorType

    zipfile_db = utils.get_test_data('prot/all.zip')
    sigs = sourmash_args.load_file_as_signatures(zipfile_db)

    # it's fine if this needs to change, but for now I want to make
    # sure that this is generator.
    assert isinstance(sigs, GeneratorType)

    sigs = list(sigs)
    assert len(sigs) == 7


def test_zipfile_load_file_as_signatures_traverse_yield_all():
    from types import GeneratorType

    zipfile_db = utils.get_test_data('prot/all.zip')
    sigs = sourmash_args.load_file_as_signatures(zipfile_db,
                                                 yield_all_files=True)

    # it's fine if this needs to change, but for now I want to make
    # sure that this is generator.
    assert isinstance(sigs, GeneratorType)

    sigs = list(sigs)
    assert len(sigs) == 8


@utils.in_tempdir
def test_zipfile_load_database_fail_if_not_zip(c):
    # fail _load_database if not .zip
    zipfile_db = utils.get_test_data('prot/all.zip')
    badname = c.output('xyz.nada')
    shutil.copyfile(zipfile_db, badname)

    with pytest.raises(ValueError) as exc:
        sigs = sourmash_args.load_file_as_signatures(badname)

    assert 'Error while reading signatures from' in str(exc.value)


def test_multi_index_search():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx1 = LinearIndex.load(sig2)
    lidx2 = LinearIndex.load(sig47)
    lidx3 = LinearIndex.load(sig63)

    # create MultiIindex with source location override
    lidx = MultiIndex([lidx1, lidx2, lidx3], ['A', None, 'C'])
    lidx = lidx.select(ksize=31)

    # now, search for sig2
    sr = lidx.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2
    assert sr[0][2] == 'A'      # source override

    # search for sig47 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss47, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss47
    assert sr[0][2] == sig47    # source was set to None, so no override
    assert sr[1][1] == ss63
    assert sr[1][2] == 'C'      # source override

    # search for sig63 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss63, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63
    assert sr[0][2] == 'C'      # source override
    assert sr[1][1] == ss47
    assert sr[1][2] == sig47    # source was set to None, so no override

    # search for sig63 with high threshold => 1 match
    sr = lidx.search(ss63, threshold=0.8)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63
    assert sr[0][2] == 'C'      # source override


def test_multi_index_gather():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx1 = LinearIndex.load(sig2)
    lidx2 = LinearIndex.load(sig47)
    lidx3 = LinearIndex.load(sig63)

    # create MultiIindex with source location override
    lidx = MultiIndex([lidx1, lidx2, lidx3], ['A', None, 'C'])
    lidx = lidx.select(ksize=31)

    matches = lidx.gather(ss2)
    assert len(matches) == 1
    assert matches[0][0] == 1.0
    assert matches[0][2] == 'A'

    matches = lidx.gather(ss47)
    assert len(matches) == 2
    assert matches[0][0] == 1.0
    assert matches[0][1] == ss47
    assert matches[0][2] == sig47     # no source override
    assert round(matches[1][0], 2) == 0.49
    assert matches[1][1] == ss63
    assert matches[1][2] == 'C'       # source override


def test_multi_index_signatures():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx1 = LinearIndex.load(sig2)
    lidx2 = LinearIndex.load(sig47)
    lidx3 = LinearIndex.load(sig63)

    # create MultiIindex with source location override
    lidx = MultiIndex([lidx1, lidx2, lidx3], ['A', None, 'C'])
    lidx = lidx.select(ksize=31)

    siglist = list(lidx.signatures())
    assert len(siglist) == 3
    assert ss2 in siglist
    assert ss47 in siglist
    assert ss63 in siglist


def test_multi_index_load_from_path():
    dirname = utils.get_test_data('prot/protein')
    mi = MultiIndex.load_from_path(dirname, force=False)

    sigs = list(mi.signatures())
    assert len(sigs) == 2


def test_multi_index_load_from_path_2():
    # only load .sig files, currently; not the databases under that directory.
    dirname = utils.get_test_data('prot')
    mi = MultiIndex.load_from_path(dirname, force=False)

    print(mi.index_list)
    print(mi.source_list)

    sigs = list(mi.signatures())
    assert len(sigs) == 7


@utils.in_tempdir
def test_multi_index_load_from_path_3(c):
    # check that force works ok on a directory
    dirname = utils.get_test_data('prot')

    count = 0
    for root, dirs, files in os.walk(dirname):
        for name in files:
            print(f"at {name}")
            fullname = os.path.join(root, name)
            copyto = c.output(f"file{count}.sig")
            shutil.copyfile(fullname, copyto)
            count += 1

    with pytest.raises(sourmash.exceptions.SourmashError):
        mi = MultiIndex.load_from_path(c.location, force=False)


@utils.in_tempdir
def test_multi_index_load_from_path_3_yield_all_true(c):
    # check that force works ok on a directory w/force=True
    dirname = utils.get_test_data('prot')

    count = 0
    for root, dirs, files in os.walk(dirname):
        for name in files:
            print(f"at {name}")
            fullname = os.path.join(root, name)
            copyto = c.output(f"file{count}.something")
            shutil.copyfile(fullname, copyto)
            count += 1

    mi = MultiIndex.load_from_path(c.location, force=True)

    print(mi.index_list)
    print(mi.source_list)

    sigs = list(mi.signatures())
    assert len(sigs) == 8


@utils.in_tempdir
def test_multi_index_load_from_path_3_yield_all_true_subdir(c):
    # check that force works ok on subdirectories
    dirname = utils.get_test_data('prot')

    target_dir = c.output("some_subdir")
    os.mkdir(target_dir)

    count = 0
    for root, dirs, files in os.walk(dirname):
        for name in files:
            print(f"at {name}")
            fullname = os.path.join(root, name)
            copyto = os.path.join(target_dir, f"file{count}.something")
            shutil.copyfile(fullname, copyto)
            count += 1

    mi = MultiIndex.load_from_path(c.location, force=True)

    print(mi.index_list)
    print(mi.source_list)

    sigs = list(mi.signatures())
    assert len(sigs) == 8


@utils.in_tempdir
def test_multi_index_load_from_path_3_sig_gz(c):
    # check that we find .sig.gz files, too
    dirname = utils.get_test_data('prot')

    count = 0
    for root, dirs, files in os.walk(dirname):
        for name in files:
            if not name.endswith('.sig'): # skip non .sig things
                continue
            print(f"at {name}")
            fullname = os.path.join(root, name)
            copyto = c.output(f"file{count}.sig.gz")
            shutil.copyfile(fullname, copyto)
            count += 1

    mi = MultiIndex.load_from_path(c.location, force=False)

    print(mi.index_list)
    print(mi.source_list)

    sigs = list(mi.signatures())
    assert len(sigs) == 6


@utils.in_tempdir
def test_multi_index_load_from_path_3_check_traverse_fn(c):
    # test the actual traverse function... eventually this test can be
    # removed, probably, as we consolidate functionality and test MultiIndex
    # better.
    dirname = utils.get_test_data('prot')
    files = list(sourmash_args.traverse_find_sigs([dirname]))
    assert len(files) == 7, files

    files = list(sourmash_args.traverse_find_sigs([dirname], True))
    assert len(files) == 20, files


def test_multi_index_load_from_path_no_exist():
    dirname = utils.get_test_data('does-not-exist')
    with pytest.raises(ValueError):
        mi = MultiIndex.load_from_path(dirname, force=True)


def test_multi_index_load_from_pathlist_no_exist():
    dirname = utils.get_test_data('does-not-exist')
    with pytest.raises(ValueError):
        mi = MultiIndex.load_from_pathlist(dirname)


@utils.in_tempdir
def test_multi_index_load_from_pathlist_1(c):
    dirname = utils.get_test_data('prot')
    files = list(sourmash_args.traverse_find_sigs([dirname]))
    assert len(files) == 7, files

    file_list = c.output('filelist.txt')

    with open(file_list, 'wt') as fp:
        print("\n".join(files), file=fp)
    mi = MultiIndex.load_from_pathlist(file_list)

    sigs = list(mi.signatures())
    assert len(sigs) == 7


@utils.in_tempdir
def test_multi_index_load_from_pathlist_2(c):
    dirname = utils.get_test_data('prot')
    files = list(sourmash_args.traverse_find_sigs([dirname], True))
    assert len(files) == 20, files

    file_list = c.output('filelist.txt')

    with open(file_list, 'wt') as fp:
        print("\n".join(files), file=fp)

    with pytest.raises(ValueError):
        mi = MultiIndex.load_from_pathlist(file_list)


@utils.in_tempdir
def test_multi_index_load_from_pathlist_3_zipfile(c):
    # can we load zipfiles in a pathlist? yes please.
    zipfile = utils.get_test_data('prot/all.zip')

    file_list = c.output('filelist.txt')

    with open(file_list, 'wt') as fp:
        print(zipfile, file=fp)

    mi = MultiIndex.load_from_pathlist(file_list)
    assert len(mi) == 7

##
## test a slightly outre version of JaccardSearch - this is a test of the
## JaccardSearch 'collect' protocol, in particular...
##

class JaccardSearchBestOnly_ButIgnore(JaccardSearch):
    "A class that ignores certain results, but still does all the pruning."
    def __init__(self, ignore_list):
        super().__init__(SearchType.JACCARD, threshold=0.1)
        self.ignore_list = ignore_list

    # a collect function that _ignores_ things in the ignore_list
    def collect(self, score, match):
        print('in collect; current threshold:', self.threshold)
        for q in self.ignore_list:
            print('ZZZ', match, match.similarity(q))
            if match.similarity(q) == 1.0:
                print('yes, found.')
                return False

        # update threshold if not perfect match, which could help prune.
        self.threshold = score
        return True


def test_linear_index_gather_ignore():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63, ksize=31)

    # construct an index...
    lidx = LinearIndex([ss2, ss47, ss63])

    # ...now search with something that should ignore sig47, the exact match.
    search_fn = JaccardSearchBestOnly_ButIgnore([ss47])

    results = list(lidx.find(search_fn, ss47))
    results = [ ss for (ss, score) in results ]

    def is_found(ss, xx):
        for q in xx:
            print(ss, ss.similarity(q))
            if ss.similarity(q) == 1.0:
                return True
        return False

    assert not is_found(ss47, results)
    assert not is_found(ss2, results)
    assert is_found(ss63, results)


def test_lca_index_gather_ignore():
    from sourmash.lca import LCA_Database

    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63, ksize=31)

    # construct an index...
    db = LCA_Database(ksize=31, scaled=1000)
    db.insert(ss2)
    db.insert(ss47)
    db.insert(ss63)

    # ...now search with something that should ignore sig47, the exact match.
    search_fn = JaccardSearchBestOnly_ButIgnore([ss47])

    results = list(db.find(search_fn, ss47))
    results = [ ss for (ss, score) in results ]

    def is_found(ss, xx):
        for q in xx:
            print(ss, ss.similarity(q))
            if ss.similarity(q) == 1.0:
                return True
        return False

    assert not is_found(ss47, results)
    assert not is_found(ss2, results)
    assert is_found(ss63, results)


def test_sbt_index_gather_ignore():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63, ksize=31)

    # construct an index...
    factory = GraphFactory(5, 100, 3)
    db = SBT(factory, d=2)

    db.insert(ss2)
    db.insert(ss47)
    db.insert(ss63)

    # ...now search with something that should ignore sig47, the exact match.
    print(f'\n** trying to ignore {ss47}')
    search_fn = JaccardSearchBestOnly_ButIgnore([ss47])

    results = list(db.find(search_fn, ss47))
    results = [ ss for (ss, score) in results ]

    def is_found(ss, xx):
        for q in xx:
            print('is found?', ss, ss.similarity(q))
            if ss.similarity(q) == 1.0:
                return True
        return False

    assert not is_found(ss47, results)
    assert not is_found(ss2, results)
    assert is_found(ss63, results)

def test_revindex_index_search():
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = RevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    # now, search for sig2
    sr = lidx.search(ss2, threshold=1.0)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2

    # search for sig47 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss47, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss47
    assert sr[1][1] == ss63

    # search for sig63 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss63, threshold=0.1)
    print([s[1].name for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63
    assert sr[1][1] == ss47

    # search for sig63 with high threshold => 1 match
    sr = lidx.search(ss63, threshold=0.8)
    print([s[1].name for s in sr])
    assert len(sr) == 1
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss63


def test_revindex_gather():
    sig2 = utils.get_test_data("2.fa.sig")
    sig47 = utils.get_test_data("47.fa.sig")
    sig63 = utils.get_test_data("63.fa.sig")

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)

    lidx = RevIndex(template=ss2.minhash)
    lidx.insert(ss2)
    lidx.insert(ss47)
    lidx.insert(ss63)

    matches = lidx.gather(ss2)
    assert len(matches) == 1
    assert matches[0][0] == 1.0
    assert matches[0][1] == ss2

    matches = lidx.gather(ss47)
    assert len(matches) == 1
    assert matches[0][0] == 1.0
    assert matches[0][1] == ss47


def test_revindex_gather_ignore():
    sig2 = utils.get_test_data('2.fa.sig')
    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss47 = sourmash.load_one_signature(sig47, ksize=31)
    ss63 = sourmash.load_one_signature(sig63, ksize=31)

    # construct an index...
    lidx = RevIndex(template=ss2.minhash, signatures=[ss2, ss47, ss63])

    # ...now search with something that should ignore sig47, the exact match.
    search_fn = JaccardSearchBestOnly_ButIgnore([ss47])

    results = list(lidx.find(search_fn, ss47))
    results = [ ss for (ss, score) in results ]

    def is_found(ss, xx):
        for q in xx:
            print(ss, ss.similarity(q))
            if ss.similarity(q) == 1.0:
                return True
        return False

    assert not is_found(ss47, results)
    assert not is_found(ss2, results)
    assert is_found(ss63, results)
