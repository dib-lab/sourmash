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
from sourmash.index import LinearIndex, MultiIndex
from sourmash.sbt import SBT, GraphFactory, Leaf
from sourmash import sourmash_args

import sourmash_tst_utils as utils


def test_simple_index(n_children):
    factory = GraphFactory(5, 100, 3)
    root = SBT(factory, d=n_children)

    leaf1 = Leaf("a", factory())
    leaf1.data.count("AAAAA")
    leaf1.data.count("AAAAT")
    leaf1.data.count("AAAAC")

    leaf2 = Leaf("b", factory())
    leaf2.data.count("AAAAA")
    leaf2.data.count("AAAAT")
    leaf2.data.count("AAAAG")

    leaf3 = Leaf("c", factory())
    leaf3.data.count("AAAAA")
    leaf3.data.count("AAAAT")
    leaf3.data.count("CAAAA")

    leaf4 = Leaf("d", factory())
    leaf4.data.count("AAAAA")
    leaf4.data.count("CAAAA")
    leaf4.data.count("GAAAA")

    leaf5 = Leaf("e", factory())
    leaf5.data.count("AAAAA")
    leaf5.data.count("AAAAT")
    leaf5.data.count("GAAAA")

    root.add_node(leaf1)
    root.add_node(leaf2)
    root.add_node(leaf3)
    root.add_node(leaf4)
    root.add_node(leaf5)

    def search_kmer(obj, seq):
        return obj.data.get(seq)

    kmers = ["AAAAA", "AAAAT", "AAAAG", "CAAAA", "GAAAA"]

    linear = LinearIndex()
    linear.insert(leaf1)
    linear.insert(leaf2)
    linear.insert(leaf3)
    linear.insert(leaf4)
    linear.insert(leaf5)

    for kmer in kmers:
        assert set(root.find(search_kmer, kmer)) == set(linear.find(search_kmer, kmer))

    print("-----")
    print([x.metadata for x in root.find(search_kmer, "AAAAA")])
    print([x.metadata for x in root.find(search_kmer, "AAAAT")])
    print([x.metadata for x in root.find(search_kmer, "AAAAG")])
    print([x.metadata for x in root.find(search_kmer, "CAAAA")])
    print([x.metadata for x in root.find(search_kmer, "GAAAA")])


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
    assert linear.filename == filename


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
def test_index_same_md5sum_zipstorage(c):
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
    assert len(sigs) == 6


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
    assert len(sigs) == 6


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
    assert len(sigs) == 6


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
    assert len(files) == 6, files

    files = list(sourmash_args.traverse_find_sigs([dirname], True))
    assert len(files) == 14, files


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
    assert len(files) == 6, files

    file_list = c.output('filelist.txt')

    with open(file_list, 'wt') as fp:
        print("\n".join(files), file=fp)
    mi = MultiIndex.load_from_pathlist(file_list)

    sigs = list(mi.signatures())
    assert len(sigs) == 6


@utils.in_tempdir
def test_multi_index_load_from_pathlist_2(c):
    dirname = utils.get_test_data('prot')
    files = list(sourmash_args.traverse_find_sigs([dirname], True))
    assert len(files) == 14, files

    file_list = c.output('filelist.txt')

    with open(file_list, 'wt') as fp:
        print("\n".join(files), file=fp)

    with pytest.raises(ValueError):
        mi = MultiIndex.load_from_pathlist(file_list)
