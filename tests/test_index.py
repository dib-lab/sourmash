"""
Tests for Index classes and subclasses.
"""
import pytest
import glob
import os
import zipfile
import shutil
import copy

import sourmash
from sourmash import load_one_signature, SourmashSignature
from sourmash.index import (LinearIndex, MultiIndex, ZipFileLinearIndex,
                            make_jaccard_search_query, CounterGather,
                            LazyLinearIndex)
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


def test_linear_index_prefetch():
    # prefetch does basic things right:
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

    # search for ss2
    results = []
    for result in lidx.prefetch(ss2, threshold_bp=0):
        results.append(result)

    assert len(results) == 1
    assert results[0].signature == ss2

    # search for ss47 - expect two results
    results = []
    for result in lidx.prefetch(ss47, threshold_bp=0):
        results.append(result)

    assert len(results) == 2
    assert results[0].signature == ss47
    assert results[1].signature == ss63


def test_linear_index_prefetch_empty():
    # check that an exception is raised upon for an empty database
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)

    lidx = LinearIndex()

    # since this is a generator, we need to actually ask for a value to
    # get exception raised.
    g = lidx.prefetch(ss2, threshold_bp=0)
    with pytest.raises(ValueError) as e:
        next(g)

    assert "no signatures to search" in str(e.value)


def test_linear_index_prefetch_lazy():
    # make sure that prefetch doesn't touch values 'til requested.
    class FakeSignature:
        @property
        def minhash(self):
            raise Exception("don't touch me!")

    sig47 = utils.get_test_data('47.fa.sig')
    sig63 = utils.get_test_data('63.fa.sig')

    ss47 = sourmash.load_one_signature(sig47)
    ss63 = sourmash.load_one_signature(sig63)
    fake = FakeSignature()

    lidx = LinearIndex()
    lidx.insert(ss47)
    lidx.insert(ss63)
    lidx.insert(fake)

    g = lidx.prefetch(ss47, threshold_bp=0)

    # first value:
    sr = next(g)
    assert sr.signature == ss47

    # second value:
    sr = next(g)
    assert sr.signature == ss63

    # third value: raises exception!
    with pytest.raises(Exception) as e:
        next(g)

    assert "don't touch me!" in str(e.value)


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
    with pytest.raises(ValueError):
        linear.gather(SourmashSignature(new_mh))

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
    with pytest.raises(ValueError):
        linear.gather(SourmashSignature(new_mh), threshold_bp=5000)

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
    with pytest.raises(ValueError):
        linear.gather(SourmashSignature(new_mh), threshold_bp=5000)


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

    print(c.last_result.out)
    print(c.last_result.err)

    assert 'no compatible signatures found in ' in c.last_result.err


def test_zipfile_API_signatures():
    # return all of the .sig and .sig.gz files in all.zip
    zipfile_db = utils.get_test_data('prot/all.zip')

    zipidx = ZipFileLinearIndex.load(zipfile_db)
    siglist = list(zipidx.signatures())
    assert len(siglist) == 7
    assert len(zipidx) == 7


def test_zipfile_bool():
    # make sure that zipfile __bool__ doesn't traverse all the signatures
    # by relying on __len__!

    # create fake class that overrides everything useful except for bool -
    class FakeZipFileLinearIndex(ZipFileLinearIndex):
        def __init__(self):
            pass

        def signatures(self):
            yield 'a'
            raise Exception("don't touch me!")

        def __len__(self):
            raise Exception("don't call len!")

    # 'bool' should not touch __len__ or a second signature
    zf = FakeZipFileLinearIndex()
    assert bool(zf)

    # __len__ should raise an exception
    with pytest.raises(Exception) as exc:
        len(zf)
    assert "don't call len!" in str(exc.value)


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
    assert len(matches) == 1
    assert matches[0][0] == 1.0
    assert matches[0][1] == ss47
    assert matches[0][2] == sig47     # no source override


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
    results = [ sr.signature for sr in results ]

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
    results = [ sr.signature for sr in results ]

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
    results = [ sr.signature for sr in results ]

    def is_found(ss, xx):
        for q in xx:
            print('is found?', ss, ss.similarity(q))
            if ss.similarity(q) == 1.0:
                return True
        return False

    assert not is_found(ss47, results)
    assert not is_found(ss2, results)
    assert is_found(ss63, results)

###
### CounterGather tests
###


def _consume_all(query_mh, counter, threshold_bp=0):
    results = []

    last_intersect_size = None
    while 1:
        result = counter.peek(query_mh, threshold_bp)
        if not result:
            break

        sr, intersect_mh = result
        print(sr.signature.name, len(intersect_mh))
        if last_intersect_size:
            assert len(intersect_mh) <= last_intersect_size

        last_intersect_size = len(intersect_mh)

        counter.consume(intersect_mh)
        query_mh.remove_many(intersect_mh.hashes)

        results.append((sr, len(intersect_mh)))

    return results


def test_counter_gather_1():
    # check a contrived set of non-overlapping gather results,
    # generated via CounterGather
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear()
    match_mh_2.add_many(range(10, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear()
    match_mh_3.add_many(range(15, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    results = _consume_all(query_ss.minhash, counter)

    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_1_b():
    # check a contrived set of somewhat-overlapping gather results,
    # generated via CounterGather. Here the overlaps are structured
    # so that the gather results are the same as those in
    # test_counter_gather_1(), even though the overlaps themselves are
    # larger.
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear()
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear()
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    results = _consume_all(query_ss.minhash, counter)

    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_1_c_with_threshold():
    # check a contrived set of somewhat-overlapping gather results,
    # generated via CounterGather. Here the overlaps are structured
    # so that the gather results are the same as those in
    # test_counter_gather_1(), even though the overlaps themselves are
    # larger.
    # use a threshold, here.

    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear()
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear()
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    results = _consume_all(query_ss.minhash, counter,
                           threshold_bp=3)

    expected = (['match1', 10],
                ['match2', 5])
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_1_d_diff_scaled():
    # test as above, but with different scaled.
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear().downsample(scaled=10)
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear().downsample(scaled=20)
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear().downsample(scaled=30)
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    results = _consume_all(query_ss.minhash, counter)

    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_1_d_diff_scaled_query():
    # test as above, but with different scaled for QUERY.
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))

    match_mh_1 = query_mh.copy_and_clear().downsample(scaled=10)
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear().downsample(scaled=20)
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear().downsample(scaled=30)
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # downsample query now -
    query_ss = SourmashSignature(query_mh.downsample(scaled=100), name='query')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    results = _consume_all(query_ss.minhash, counter)

    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_1_e_abund_query():
    # test as above, but abund query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1, track_abundance=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear().flatten()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear().flatten()
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear().flatten()
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    # must flatten before peek!
    results = _consume_all(query_ss.minhash.flatten(), counter)

    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_1_f_abund_match():
    # test as above, but abund query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1, track_abundance=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh.flatten(), name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear()
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear()
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(match_ss_1)
    counter.add(match_ss_2)
    counter.add(match_ss_3)

    # must flatten before peek!
    results = _consume_all(query_ss.minhash.flatten(), counter)

    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_2():
    # check basic set of gather results on semi-real data,
    # generated via CounterGather
    testdata_combined = utils.get_test_data('gather/combined.sig')
    testdata_glob = utils.get_test_data('gather/GCF*.sig')
    testdata_sigs = glob.glob(testdata_glob)

    query_ss = sourmash.load_one_signature(testdata_combined, ksize=21)
    subject_sigs = [ (sourmash.load_one_signature(t, ksize=21), t)
                     for t in testdata_sigs ]

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    for ss, loc in subject_sigs:
        counter.add(ss, loc)

    results = _consume_all(query_ss.minhash, counter)

    expected = (['NC_003198.1', 487],
                ['NC_000853.1', 192],
                ['NC_011978.1', 169],
                ['NC_002163.1', 157],
                ['NC_003197.2', 152],
                ['NC_009486.1', 92],
                ['NC_006905.1', 76],
                ['NC_011080.1', 59],
                ['NC_011274.1', 42],
                ['NC_006511.1', 31],
                ['NC_011294.1', 7],
                ['NC_004631.1', 2])
    assert len(results) == len(expected)

    for (sr, size), (exp_name, exp_size) in zip(results, expected):
        sr_name = sr.signature.name.split()[0]
        print(sr_name, size)

        assert sr_name == exp_name
        assert size == exp_size


def test_counter_gather_exact_match():
    # query == match
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(query_ss, 'somewhere over the rainbow')

    results = _consume_all(query_ss.minhash, counter)
    assert len(results) == 1
    (sr, intersect_mh) = results[0]

    assert sr.score == 1.0
    assert sr.signature == query_ss
    assert sr.location == 'somewhere over the rainbow'


def test_counter_gather_add_after_peek():
    # cannot add after peek or consume
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(query_ss, 'somewhere over the rainbow')

    counter.peek(query_ss.minhash)

    with pytest.raises(ValueError):
        counter.add(query_ss, "try again")


def test_counter_gather_add_after_consume():
    # cannot add after peek or consume
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(query_ss, 'somewhere over the rainbow')

    counter.consume(query_ss.minhash)

    with pytest.raises(ValueError):
        counter.add(query_ss, "try again")


def test_counter_gather_consume_empty_intersect():
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(query_ss, 'somewhere over the rainbow')

    # nothing really happens here :laugh:, just making sure there's no error
    counter.consume(query_ss.minhash.copy_and_clear())


def test_counter_gather_empty_initial_query():
    # check empty initial query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(match_ss_1, require_overlap=False)

    assert counter.peek(query_ss.minhash) == []


def test_counter_gather_num_query():
    # check num query
    query_mh = sourmash.MinHash(n=500, ksize=31)
    query_mh.add_many(range(0, 10))
    query_ss = SourmashSignature(query_mh, name='query')

    with pytest.raises(ValueError):
        counter = CounterGather(query_ss.minhash)


def test_counter_gather_empty_cur_query():
    # test empty cur query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(query_ss, 'somewhere over the rainbow')

    cur_query_mh = query_ss.minhash.copy_and_clear()
    results = _consume_all(cur_query_mh, counter)
    assert results == []


def test_counter_gather_add_num_matchy():
    # test add num query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh = sourmash.MinHash(n=500, ksize=31)
    match_mh.add_many(range(0, 20))
    match_ss = SourmashSignature(match_mh, name='query')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    with pytest.raises(ValueError):
        counter.add(match_ss, 'somewhere over the rainbow')


def test_counter_gather_bad_cur_query():
    # test cur query that is not subset of original query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(query_ss, 'somewhere over the rainbow')

    cur_query_mh = query_ss.minhash.copy_and_clear()
    cur_query_mh.add_many(range(20, 30))
    with pytest.raises(ValueError):
        counter.peek(cur_query_mh)


def test_counter_gather_add_no_overlap():
    # check adding match with no overlap w/query
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 10))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(10, 20))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    with pytest.raises(ValueError):
        counter.add(match_ss_1)

    assert counter.peek(query_ss.minhash) == []


def test_counter_gather_big_threshold():
    # check 'peek' with a huge threshold
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(match_ss_1)

    # impossible threshold:
    threshold_bp=30*query_ss.minhash.scaled
    results = counter.peek(query_ss.minhash, threshold_bp=threshold_bp)
    assert results == []


def test_counter_gather_empty_counter():
    # check empty counter
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_ss = SourmashSignature(query_mh, name='query')

    # empty counter!
    counter = CounterGather(query_ss.minhash)

    assert counter.peek(query_ss.minhash) == []


def test_counter_gather_3_test_consume():
    # open-box testing of consume(...)
    query_mh = sourmash.MinHash(n=0, ksize=31, scaled=1)
    query_mh.add_many(range(0, 20))
    query_ss = SourmashSignature(query_mh, name='query')

    match_mh_1 = query_mh.copy_and_clear()
    match_mh_1.add_many(range(0, 10))
    match_ss_1 = SourmashSignature(match_mh_1, name='match1')

    match_mh_2 = query_mh.copy_and_clear()
    match_mh_2.add_many(range(7, 15))
    match_ss_2 = SourmashSignature(match_mh_2, name='match2')

    match_mh_3 = query_mh.copy_and_clear()
    match_mh_3.add_many(range(13, 17))
    match_ss_3 = SourmashSignature(match_mh_3, name='match3')

    # load up the counter
    counter = CounterGather(query_ss.minhash)
    counter.add(match_ss_1, 'loc a')
    counter.add(match_ss_2, 'loc b')
    counter.add(match_ss_3, 'loc c')

    ### ok, dig into actual counts...
    import pprint
    pprint.pprint(counter.counter)
    pprint.pprint(counter.siglist)
    pprint.pprint(counter.locations)

    assert counter.siglist == [ match_ss_1, match_ss_2, match_ss_3 ]
    assert counter.locations == ['loc a', 'loc b', 'loc c']
    assert list(counter.counter.items()) == [(0, 10), (1, 8), (2, 4)]

    ## round 1

    cur_query = copy.copy(query_ss.minhash)
    (sr, intersect_mh) = counter.peek(cur_query)
    assert sr.signature == match_ss_1
    assert len(intersect_mh) == 10
    assert cur_query == query_ss.minhash

    counter.consume(intersect_mh)
    assert counter.siglist == [ match_ss_1, match_ss_2, match_ss_3 ]
    assert counter.locations == ['loc a', 'loc b', 'loc c']
    assert list(counter.counter.items()) == [(1, 5), (2, 4)]

    ### round 2

    cur_query.remove_many(intersect_mh.hashes)
    (sr, intersect_mh) = counter.peek(cur_query)
    assert sr.signature == match_ss_2
    assert len(intersect_mh) == 5
    assert cur_query != query_ss.minhash

    counter.consume(intersect_mh)
    assert counter.siglist == [ match_ss_1, match_ss_2, match_ss_3 ]
    assert counter.locations == ['loc a', 'loc b', 'loc c']
    assert list(counter.counter.items()) == [(2, 2)]

    ## round 3

    cur_query.remove_many(intersect_mh.hashes)
    (sr, intersect_mh) = counter.peek(cur_query)
    assert sr.signature == match_ss_3
    assert len(intersect_mh) == 2
    assert cur_query != query_ss.minhash

    counter.consume(intersect_mh)
    assert counter.siglist == [ match_ss_1, match_ss_2, match_ss_3 ]
    assert counter.locations == ['loc a', 'loc b', 'loc c']
    assert list(counter.counter.items()) == []

    ## round 4 - nothing left!

    cur_query.remove_many(intersect_mh.hashes)
    results = counter.peek(cur_query)
    assert not results

    counter.consume(intersect_mh)
    assert counter.siglist == [ match_ss_1, match_ss_2, match_ss_3 ]
    assert counter.locations == ['loc a', 'loc b', 'loc c']
    assert list(counter.counter.items()) == []


def test_lazy_index_1():
    # test some basic features of LazyLinearIndex
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

    lazy = LazyLinearIndex(lidx)
    lazy2 = lazy.select(ksize=31)
    assert len(list(lazy2.signatures())) == 3

    results = lazy2.search(ss2, threshold=1.0)
    assert len(results) == 1
    assert results[0].signature == ss2


def test_lazy_index_2():
    # test laziness by adding a signature that raises an exception when
    # touched.

    class FakeSignature:
        @property
        def minhash(self):
            raise Exception("don't touch me!")

    lidx = LinearIndex()
    lidx.insert(FakeSignature())

    lazy = LazyLinearIndex(lidx)
    lazy2 = lazy.select(ksize=31)

    sig_iter = lazy2.signatures()
    with pytest.raises(Exception) as e:
        list(sig_iter)

    assert str(e.value) == "don't touch me!"


def test_lazy_index_3():
    # make sure that you can't do multiple _incompatible_ selects.
    class FakeSignature:
        @property
        def minhash(self):
            raise Exception("don't touch me!")

    lidx = LinearIndex()
    lidx.insert(FakeSignature())

    lazy = LazyLinearIndex(lidx)
    lazy2 = lazy.select(ksize=31)
    with pytest.raises(ValueError) as e:
        lazy3 = lazy2.select(ksize=21)

    assert str(e.value) == "cannot select on two different values for ksize"


def test_lazy_index_4_bool():
    # test some basic features of LazyLinearIndex
    sig2 = utils.get_test_data('2.fa.sig')
    ss2 = sourmash.load_one_signature(sig2, ksize=31)

    # test bool false/true
    lidx = LinearIndex()
    lazy = LazyLinearIndex(lidx)
    assert not lazy

    lidx.insert(ss2)
    assert lazy


def test_lazy_index_5_len():
    # test some basic features of LazyLinearIndex
    lidx = LinearIndex()
    lazy = LazyLinearIndex(lidx)

    with pytest.raises(NotImplementedError):
        len(lazy)


def test_lazy_index_wraps_multiindex_location():
    sigdir = utils.get_test_data('prot/protein/')
    sigzip = utils.get_test_data('prot/protein.zip')
    siglca = utils.get_test_data('prot/protein.lca.json.gz')
    sigsbt = utils.get_test_data('prot/protein.sbt.zip')

    db_paths = (sigdir, sigzip, siglca, sigsbt)
    dbs = [ sourmash.load_file_as_index(db_path) for db_path in db_paths ]

    mi = MultiIndex(dbs, db_paths)
    lazy = LazyLinearIndex(mi)

    mi2 = mi.select(moltype='protein')
    lazy2 = lazy.select(moltype='protein')

    for (ss_tup, ss_lazy_tup) in zip(mi2.signatures_with_location(),
                                     lazy2.signatures_with_location()):
        assert ss_tup == ss_lazy_tup
