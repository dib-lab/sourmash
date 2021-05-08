"Test SBT code."
import json
import shutil
import os

import pytest

import sourmash
from sourmash import (load_one_signature, SourmashSignature,
                      load_file_as_signatures)
from sourmash.exceptions import IndexNotSupported
from sourmash.sbt import SBT, GraphFactory, Leaf, Node
from sourmash.sbtmh import (SigLeaf, load_sbt_index)
from sourmash.sbt_storage import (FSStorage, RedisStorage,
                                  IPFSStorage, ZipStorage)
from sourmash.search import make_jaccard_search_query

import sourmash_tst_utils as utils


def test_simple(n_children):
    factory = GraphFactory(5, 100, 3)
    root = SBT(factory, d=n_children)

    leaf1 = Leaf("a", factory())
    leaf1.data.count('AAAAA')
    leaf1.data.count('AAAAT')
    leaf1.data.count('AAAAC')

    leaf2 = Leaf("b", factory())
    leaf2.data.count('AAAAA')
    leaf2.data.count('AAAAT')
    leaf2.data.count('AAAAG')

    leaf3 = Leaf("c", factory())
    leaf3.data.count('AAAAA')
    leaf3.data.count('AAAAT')
    leaf3.data.count('CAAAA')

    leaf4 = Leaf("d", factory())
    leaf4.data.count('AAAAA')
    leaf4.data.count('CAAAA')
    leaf4.data.count('GAAAA')

    leaf5 = Leaf("e", factory())
    leaf5.data.count('AAAAA')
    leaf5.data.count('AAAAT')
    leaf5.data.count('GAAAA')

    root.add_node(leaf1)
    root.add_node(leaf2)
    root.add_node(leaf3)
    root.add_node(leaf4)
    root.add_node(leaf5)

    # return True if leaf node contains nodegraph w/kmer
    def search_kmer(leaf, kmer):
        return leaf.data.get(kmer)

    leaves = [leaf1, leaf2, leaf3, leaf4, leaf5 ]
    kmers = [ "AAAAA", "AAAAT", "AAAAG", "CAAAA", "GAAAA" ]

    # define an exhaustive search function that looks in all the leaf nodes.
    def search_kmer_in_list(kmer):
        x = []
        for l in leaves:
            if l.data.get(kmer):
                x.append(l)

        return set(x)

    # for all k-mers, ensure that tree._find_nodes matches the exhaustive
    # search.
    for kmer in kmers:
        assert set(root._find_nodes(search_kmer, kmer)) == search_kmer_in_list(kmer)

    print('-----')
    print([ x.metadata for x in root._find_nodes(search_kmer, "AAAAA") ])
    print([ x.metadata for x in root._find_nodes(search_kmer, "AAAAT") ])
    print([ x.metadata for x in root._find_nodes(search_kmer, "AAAAG") ])
    print([ x.metadata for x in root._find_nodes(search_kmer, "CAAAA") ])
    print([ x.metadata for x in root._find_nodes(search_kmer, "GAAAA") ])

    # save SBT to a directory and then reload
    with utils.TempDirectory() as location:
        root.save(os.path.join(location, 'demo'))
        root = SBT.load(os.path.join(location, 'demo'))

        for kmer in kmers:
            new_result = {str(r) for r in root._find_nodes(search_kmer, kmer)}
            print(*new_result, sep='\n')

            assert new_result == {str(r) for r in search_kmer_in_list(kmer)}


def test_longer_search(n_children):
    ksize = 5
    factory = GraphFactory(ksize, 100, 3)
    root = SBT(factory, d=n_children)

    leaf1 = Leaf("a", factory())
    leaf1.data.count('AAAAA')
    leaf1.data.count('AAAAT')
    leaf1.data.count('AAAAC')

    leaf2 = Leaf("b", factory())
    leaf2.data.count('AAAAA')
    leaf2.data.count('AAAAT')
    leaf2.data.count('AAAAG')

    leaf3 = Leaf("c", factory())
    leaf3.data.count('AAAAA')
    leaf3.data.count('AAAAT')
    leaf3.data.count('CAAAA')

    leaf4 = Leaf("d", factory())
    leaf4.data.count('AAAAA')
    leaf4.data.count('CAAAA')
    leaf4.data.count('GAAAA')

    leaf5 = Leaf("e", factory())
    leaf5.data.count('AAAAA')
    leaf5.data.count('AAAAT')
    leaf5.data.count('GAAAA')

    root.add_node(leaf1)
    root.add_node(leaf2)
    root.add_node(leaf3)
    root.add_node(leaf4)
    root.add_node(leaf5)

    def kmers(k, seq):
        for start in range(len(seq) - k + 1):
            yield seq[start:start + k]

    def search_transcript(node, seq, threshold):
        presence = [ node.data.get(kmer) for kmer in kmers(ksize, seq) ]
        if sum(presence) >= int(threshold * (len(seq) - ksize + 1)):
            return 1
        return 0

    try1 = [ x.metadata for x in root._find_nodes(search_transcript, "AAAAT", 1.0) ]
    assert set(try1) == set([ 'a', 'b', 'c', 'e' ]), try1 # no 'd'

    try2 = [ x.metadata for x in root._find_nodes(search_transcript, "GAAAAAT", 0.6) ]
    assert set(try2) == set([ 'a', 'b', 'c', 'd', 'e' ])

    try3 = [ x.metadata for x in root._find_nodes(search_transcript, "GAAAA", 1.0) ]
    assert set(try3) == set([ 'd', 'e' ]), try3


#@pytest.mark.parametrize("old_version", ["v1", "v2", "v3", "v4", "v5"])
@pytest.mark.parametrize("old_version", ["v3", "v4", "v5"])
def test_tree_old_load(old_version):
    tree_old = SBT.load(utils.get_test_data('{}.sbt.json'.format(old_version)),
                       leaf_loader=SigLeaf.load)

    tree_cur = SBT.load(utils.get_test_data('v6.sbt.json'),
                        leaf_loader=SigLeaf.load)

    testdata1 = utils.get_test_data(utils.SIG_FILES[0])
    to_search = load_one_signature(testdata1)

    print(list(tree_old.leaves()))

    # note: earlier versions of this test did containment on
    # the num MinHash in `to_search`, which doesn't work properly.
    # (See test_sbt_no_containment_on_num for test). So, to
    # fix the test for the new search API, we had to adjust
    # the threshold.
    search_obj = make_jaccard_search_query(threshold=0.05)
    results_old = {str(s.signature) for s in tree_old.find(search_obj, to_search)}
    results_cur = {str(s.signature) for s in tree_cur.find(search_obj, to_search)}

    assert results_old == results_cur
    assert len(results_old) == 4


def test_load_future(tmpdir):
    with open(str(tmpdir.join("v9999.sbt.json")), 'w') as f:
        json.dump({'version': 9999}, f)

    with pytest.raises(IndexNotSupported) as excinfo:
        SBT.load(str(tmpdir.join("v9999.sbt.json")))

    assert "index format is not supported" in str(excinfo.value)


def test_tree_save_load(n_children):
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=n_children)

    for f in utils.SIG_FILES:
        sig = load_one_signature(utils.get_test_data(f))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree.add_node(leaf)
        to_search = leaf

    print('*' * 60)
    print("{}:".format(to_search.metadata))
    search_obj = make_jaccard_search_query(threshold=0.1)
    old_result = {str(s.signature) for s in tree.find(search_obj, to_search.data)}
    print(*old_result, sep='\n')

    with utils.TempDirectory() as location:
        tree.save(os.path.join(location, 'demo'))
        tree = SBT.load(os.path.join(location, 'demo'),
                        leaf_loader=SigLeaf.load)

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        search_obj = make_jaccard_search_query(threshold=0.1)
        new_result = {str(s.signature) for s in tree.find(search_obj, to_search.data)}
        print(*new_result, sep='\n')

        assert old_result == new_result


def test_search_minhashes():
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory)

    n_leaves = 0
    for f in utils.SIG_FILES:
        sig = load_one_signature(utils.get_test_data(f))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree.add_node(leaf)

    to_search = next(iter(tree.leaves()))

    # this fails if 'search_obj' is calc containment and not similarity.
    search_obj = make_jaccard_search_query(threshold=0.08)
    results = tree.find(search_obj, to_search.data)
    for sr in results:
        assert to_search.data.jaccard(sr.signature) >= 0.08

    print(results)


def test_binary_nary_tree():
    factory = GraphFactory(31, 1e5, 4)
    trees = {}
    trees[2] = SBT(factory)
    trees[5] = SBT(factory, d=5)
    trees[10] = SBT(factory, d=10)

    n_leaves = 0
    for f in utils.SIG_FILES:
        sig = load_one_signature(utils.get_test_data(f))
        leaf = SigLeaf(os.path.basename(f), sig)
        for tree in trees.values():
            tree.add_node(leaf)
        to_search = leaf
        n_leaves += 1

    assert all([len(list(t.leaves())) == n_leaves for t in trees.values()])

    results = {}
    print('*' * 60)
    print("{}:".format(to_search.metadata))
    for d, tree in trees.items():
        search_obj = make_jaccard_search_query(threshold=0.1)
        results[d] = {str(s.signature) for s in tree.find(search_obj, to_search.data)}
    print(*results[2], sep='\n')

    assert results[2] == results[5]
    assert results[5] == results[10]


def test_sbt_combine(n_children):
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=n_children)
    tree_1 = SBT(factory, d=n_children)
    tree_2 = SBT(factory, d=n_children)

    n_leaves = 0
    for f in utils.SIG_FILES:
        sig = load_one_signature(utils.get_test_data(f))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree.add_node(leaf)
        if n_leaves < 4:
            tree_1.add_node(leaf)
        else:
            tree_2.add_node(leaf)
        n_leaves += 1

    tree_1.combine(tree_2)

    t1_leaves = {str(l) for l in tree_1.leaves()}
    t_leaves = {str(l) for l in tree.leaves()}

    assert len(t1_leaves) == n_leaves
    assert len(t_leaves) == len(t1_leaves)
    assert t1_leaves == t_leaves

    to_search = load_one_signature(utils.get_test_data(utils.SIG_FILES[0]))
    search_obj = make_jaccard_search_query(threshold=0.1)
    t1_result = {str(s.signature) for s in tree_1.find(search_obj, to_search)}
    tree_result = {str(s.signature) for s in tree.find(search_obj, to_search)}
    assert t1_result == tree_result

    # TODO: save and load both trees

    # check if adding a new node will use the next empty position
    next_empty = 0
    for n, (d, _) in enumerate(tree_1):
        if n != d:
            next_empty = n
            break
    if not next_empty:
        next_empty = n + 1

    tree_1.add_node(SigLeaf(to_search.name, to_search))
    assert tree_1.next_node == next_empty


def test_sbt_fsstorage():
    factory = GraphFactory(31, 1e5, 4)
    with utils.TempDirectory() as location:
        tree = SBT(factory)

        for f in utils.SIG_FILES:
            sig = load_one_signature(utils.get_test_data(f))

            leaf = SigLeaf(os.path.basename(f), sig)
            tree.add_node(leaf)
            to_search = leaf

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        search_obj = make_jaccard_search_query(threshold=0.1)
        old_result = {str(s.signature) for s in tree.find(search_obj, to_search.data)}
        print(*old_result, sep='\n')

        with FSStorage(location, '.fstree') as storage:
            tree.save(os.path.join(location, 'tree.sbt.json'), storage=storage)

        tree = SBT.load(os.path.join(location, 'tree.sbt.json'), leaf_loader=SigLeaf.load)
        print('*' * 60)
        print("{}:".format(to_search.metadata))
        search_obj = make_jaccard_search_query(threshold=0.1)
        new_result = {str(s.signature) for s in tree.find(search_obj, to_search.data)}
        print(*new_result, sep='\n')

        assert old_result == new_result

        assert os.path.exists(os.path.join(location, tree.storage.subdir))
        assert os.path.exists(os.path.join(location, '.fstree'))


def test_sbt_zipstorage(tmpdir):
    # create tree, save to a zip, then load and search.
    factory = GraphFactory(31, 1e5, 4)

    tree = SBT(factory)

    for f in utils.SIG_FILES:
        sig = next(load_file_as_signatures(utils.get_test_data(f)))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree.add_node(leaf)
        to_search = leaf

    print('*' * 60)
    print("{}:".format(to_search.metadata))
    search_obj = make_jaccard_search_query(threshold=0.1)
    old_result = {str(s.signature) for s in tree.find(search_obj, to_search.data)}
    print(*old_result, sep='\n')

    with ZipStorage(str(tmpdir.join("tree.sbt.zip"))) as storage:
        tree.save(str(tmpdir.join("tree.sbt.json")), storage=storage)

    with ZipStorage(str(tmpdir.join("tree.sbt.zip"))) as storage:
        tree = SBT.load(str(tmpdir.join("tree.sbt.json")),
                        leaf_loader=SigLeaf.load,
                        storage=storage)

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        search_obj = make_jaccard_search_query(threshold=0.1)
        new_result = {str(s.signature) for s in tree.find(search_obj, to_search.data)}
        print(*new_result, sep='\n')

        assert old_result == new_result


def test_sbt_ipfsstorage():
    ipfshttpclient = pytest.importorskip('ipfshttpclient')

    factory = GraphFactory(31, 1e5, 4)
    with utils.TempDirectory() as location:
        tree = SBT(factory)

        for f in utils.SIG_FILES:
            sig = load_one_signature(utils.get_test_data(f))

            leaf = SigLeaf(os.path.basename(f), sig)
            tree.add_node(leaf)
            to_search = leaf

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        search_obj = make_jaccard_search_query(threshold=0.1)
        old_result = {str(s.signature) for s in tree.find(search_obj, to_search.data)}
        print(*old_result, sep='\n')

        try:
            with IPFSStorage() as storage:
                tree.save(os.path.join(location, 'tree.sbt.json'), storage=storage)
        except ipfshttpclient.exceptions.ConnectionError:
            pytest.xfail("ipfs not installed/functioning probably")

        with IPFSStorage() as storage:
            tree = SBT.load(os.path.join(location, 'tree.sbt.json'),
                            leaf_loader=SigLeaf.load,
                            storage=storage)

            print('*' * 60)
            print("{}:".format(to_search.metadata))
            search_obj = make_jaccard_search_query(threshold=0.1)
            new_result = {str(s.signature) for s in tree.find(search_obj, to_search.data)}
            print(*new_result, sep='\n')

            assert old_result == new_result


def test_sbt_redisstorage():
    redis = pytest.importorskip('redis')
    factory = GraphFactory(31, 1e5, 4)
    with utils.TempDirectory() as location:
        tree = SBT(factory)

        for f in utils.SIG_FILES:
            sig = load_one_signature(utils.get_test_data(f))

            leaf = SigLeaf(os.path.basename(f), sig)
            tree.add_node(leaf)
            to_search = leaf

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        search_obj = make_jaccard_search_query(threshold=0.1)
        old_result = {str(s.signature) for s in tree.find(search_obj, to_search.data)}
        print(*old_result, sep='\n')

        try:
            with RedisStorage() as storage:
                tree.save(os.path.join(location, 'tree.sbt.json'), storage=storage)
        except redis.exceptions.ConnectionError:
            pytest.xfail("Couldn't connect to redis server")

        with RedisStorage() as storage:
            tree = SBT.load(os.path.join(location, 'tree.sbt.json'),
                            leaf_loader=SigLeaf.load,
                            storage=storage)

            print('*' * 60)
            print("{}:".format(to_search.metadata))
            search_obj = make_jaccard_search_query(threshold=0.1)
            new_result = {str(s.signature) for s in tree.find(search_obj, to_search.data)}
            print(*new_result, sep='\n')

            assert old_result == new_result


def test_save_zip(tmpdir):
    # load from zipped SBT, save to zipped SBT, and then search.
    testdata = utils.get_test_data("v6.sbt.zip")
    testsbt = tmpdir.join("v6.sbt.zip")
    newsbt = tmpdir.join("new.sbt.zip")

    shutil.copyfile(testdata, str(testsbt))

    tree = SBT.load(str(testsbt), leaf_loader=SigLeaf.load)
    tree.save(str(newsbt))
    assert newsbt.exists()

    new_tree = SBT.load(str(newsbt), leaf_loader=SigLeaf.load)
    assert isinstance(new_tree.storage, ZipStorage)
    assert new_tree.storage.list_sbts() == ['new.sbt.json']

    to_search = load_one_signature(utils.get_test_data(utils.SIG_FILES[0]))

    print("*" * 60)
    print("{}:".format(to_search))
    search_obj = make_jaccard_search_query(threshold=0.1)
    old_result = {str(s.signature) for s in tree.find(search_obj, to_search)}
    new_result = {str(s.signature) for s in new_tree.find(search_obj, to_search)}
    print(*new_result, sep="\n")

    assert old_result == new_result
    assert len(new_result) == 2


def test_load_zip(tmpdir):
    # search zipped SBT
    testdata = utils.get_test_data("v6.sbt.zip")
    testsbt = tmpdir.join("v6.sbt.zip")

    shutil.copyfile(testdata, str(testsbt))

    tree = SBT.load(str(testsbt), leaf_loader=SigLeaf.load)

    to_search = load_one_signature(utils.get_test_data(utils.SIG_FILES[0]))

    print("*" * 60)
    print("{}:".format(to_search))
    search_obj = make_jaccard_search_query(threshold=0.1)
    new_result = {str(s.signature) for s in tree.find(search_obj, to_search)}
    print(*new_result, sep="\n")
    assert len(new_result) == 2


def test_load_zip_uncompressed(tmpdir):
    # uncompress zipped SBT into a tmpdir and search unpacked SBT
    import zipfile

    testdata = utils.get_test_data("v6.sbt.zip")
    testsbt = tmpdir.join("v6.sbt.json")

    with zipfile.ZipFile(testdata, 'r') as z:
        z.extractall(str(tmpdir))

    tree = SBT.load(str(testsbt), leaf_loader=SigLeaf.load)

    to_search = load_one_signature(utils.get_test_data(utils.SIG_FILES[0]))

    print("*" * 60)
    print("{}:".format(to_search))
    search_obj = make_jaccard_search_query(threshold=0.1)
    new_result = {str(s.signature) for s in tree.find(search_obj, to_search)}
    print(*new_result, sep="\n")
    assert len(new_result) == 2


def test_tree_repair():
    tree_repair = SBT.load(utils.get_test_data('leaves.sbt.json'),
                           leaf_loader=SigLeaf.load)

    tree_cur = SBT.load(utils.get_test_data('v3.sbt.json'),
                        leaf_loader=SigLeaf.load)

    testdata1 = utils.get_test_data(utils.SIG_FILES[0])
    to_search = load_one_signature(testdata1)

    search_obj = make_jaccard_search_query(threshold=0.1)
    results_repair = {str(s.signature) for s in tree_repair.find(search_obj, to_search)}
    results_cur = {str(s.signature) for s in tree_cur.find(search_obj, to_search)}

    assert results_repair == results_cur
    assert len(results_repair) == 2


def test_tree_repair_insert():
    tree_repair = SBT.load(utils.get_test_data('leaves.sbt.json'),
                           leaf_loader=SigLeaf.load)

    for f in utils.SIG_FILES:
        sig = load_one_signature(utils.get_test_data(f))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree_repair.add_node(leaf)

    for pos, node in tree_repair:
        # Every parent of a node must be an internal node (and not a leaf),
        # except for node 0 (the root), whose parent is None.
        if pos != 0:
            assert isinstance(tree_repair.parent(pos).node, Node)

        # Leaf nodes can't have children
        if isinstance(node, Leaf):
            assert all(c.node is None for c in tree_repair.children(pos))


def test_save_sparseness(n_children):
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=n_children)

    for f in utils.SIG_FILES:
        sig = load_one_signature(utils.get_test_data(f))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree.add_node(leaf)
        to_search = leaf

    print('*' * 60)
    print("{}:".format(to_search.metadata))

    search_obj = make_jaccard_search_query(threshold=0.1)
    old_result = {str(s.signature) for s in tree.find(search_obj, to_search.data)}
    print(*old_result, sep='\n')

    with utils.TempDirectory() as location:
        tree.save(os.path.join(location, 'demo'), sparseness=1.0)
        tree_loaded = SBT.load(os.path.join(location, 'demo'),
                               leaf_loader=SigLeaf.load)
        assert all(not isinstance(n, Node) for _, n in tree_loaded)

        print('*' * 60)
        print("{}:".format(to_search.metadata))
        new_result = {str(s.signature) for s in tree_loaded.find(search_obj,
                                                       to_search.data)}
        print(*new_result, sep='\n')

        assert old_result == new_result

        for pos, node in tree_loaded:
            # Every parent of a node must be an internal node (and not a leaf),
            # except for node 0 (the root), whose parent is None.
            if pos != 0:
                assert isinstance(tree_loaded.parent(pos).node, Node)

            # Leaf nodes can't have children
            if isinstance(node, Leaf):
                assert all(c.node is None for c in tree_loaded.children(pos))


def test_sbt_as_index_select():
    # test 'select' method from Index base class.
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=2)

    sig47 = load_one_signature(utils.get_test_data('47.fa.sig'))
    sig63 = load_one_signature(utils.get_test_data('63.fa.sig'))

    tree.insert(sig47)
    tree.insert(sig63)

    xx = tree.select(ksize=31)
    assert xx == tree

    xx = tree.select(moltype='DNA')
    assert xx == tree

    with pytest.raises(ValueError):
        tree.select(ksize=21)

    with pytest.raises(ValueError):
        tree.select(moltype='protein')


def test_sbt_as_index_signatures():
    # test 'signatures' method from Index base class.
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=2)

    sig47 = load_one_signature(utils.get_test_data('47.fa.sig'))
    sig63 = load_one_signature(utils.get_test_data('63.fa.sig'))

    tree.insert(sig47)
    tree.insert(sig63)

    xx = list(tree.signatures())
    assert len(xx) == 2

    assert sig47 in xx
    assert sig63 in xx


def test_sbt_gather_threshold_1():
    # test gather() method, in some detail
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=2)

    sig2 = load_one_signature(utils.get_test_data('2.fa.sig'), ksize=31)
    sig47 = load_one_signature(utils.get_test_data('47.fa.sig'), ksize=31)
    sig63 = load_one_signature(utils.get_test_data('63.fa.sig'), ksize=31)

    tree.insert(sig47)
    tree.insert(sig63)
    tree.insert(sig2)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(sig2.minhash.hashes.keys()))
    new_mh = sig2.minhash.copy_and_clear()

    # query with empty hashes
    assert not new_mh
    with pytest.raises(ValueError):
        tree.gather(SourmashSignature(new_mh))

    # add one hash
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 1

    results = tree.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig == sig2
    assert name is None

    # check with a threshold -> should be no results.
    with pytest.raises(ValueError):
        tree.gather(SourmashSignature(new_mh), threshold_bp=5000)

    # add three more hashes => length of 4
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 4

    results = tree.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig == sig2
    assert name is None

    # check with a too-high threshold -> should be no results.
    print('len mh', len(new_mh))
    with pytest.raises(ValueError):
        tree.gather(SourmashSignature(new_mh), threshold_bp=5000)


def test_sbt_gather_threshold_5():
    # test gather() method above threshold
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=2)

    sig2 = load_one_signature(utils.get_test_data('2.fa.sig'), ksize=31)
    sig47 = load_one_signature(utils.get_test_data('47.fa.sig'), ksize=31)
    sig63 = load_one_signature(utils.get_test_data('63.fa.sig'), ksize=31)

    tree.insert(sig47)
    tree.insert(sig63)
    tree.insert(sig2)

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
    results = tree.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig == sig2
    assert name is None

    # now, check with a threshold_bp that should be meet-able.
    results = tree.gather(SourmashSignature(new_mh), threshold_bp=5000)
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig == sig2
    assert name is None


@utils.in_tempdir
def test_gather_single_return(c):
    # test gather() number of returns
    sig2file = utils.get_test_data('2.fa.sig')
    sig47file = utils.get_test_data('47.fa.sig')
    sig63file = utils.get_test_data('63.fa.sig')

    sig2 = load_one_signature(sig2file, ksize=31)
    sig47 = load_one_signature(sig47file, ksize=31)
    sig63 = load_one_signature(sig63file, ksize=31)

    # construct LCA Database
    factory = GraphFactory(31, 1e5, 4)
    tree = SBT(factory, d=2)

    tree.insert(sig2)
    tree.insert(sig47)
    tree.insert(sig63)

    # now, run gather. how many results do we get, and are they in the
    # right order?
    results = tree.gather(sig63)
    print(len(results))
    assert len(results) == 1
    assert results[0][0] == 1.0


@utils.in_tempdir
def test_sbt_protein_command_index(c):
    # test command-line creation of SBT database with protein sigs
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    db_out = c.output('protein.sbt.zip')

    c.run_sourmash('index', db_out, sigfile1, sigfile2,
                   '--scaled', '100', '-k', '19', '--protein')

    db2 = load_sbt_index(db_out)

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0, ignore_abundance=True,
                         do_containment=False, best_only=False)
    assert len(results) == 2

    results = db2.gather(sig2)
    assert results[0][0] == 1.0
    assert results[0][2] == db2._location
    assert results[0][2] == db_out


@utils.in_tempdir
def test_sbt_protein_search_no_threshold(c):
    # test the '.search' method on SBTs w/no threshold
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    db_out = c.output('protein.sbt.zip')

    c.run_sourmash('index', db_out, sigfile1, sigfile2,
                   '--scaled', '100', '-k', '19', '--protein')

    db2 = load_sbt_index(db_out)

    sig1 = sourmash.load_one_signature(sigfile1)

    # and search, gather
    with pytest.raises(TypeError) as exc:
        results = db2.search(sig1)
    assert "'search' requires 'threshold'" in str(exc)


@utils.in_thisdir
def test_sbt_protein_command_search(c):
    # test command-line search/gather of SBT database with protein sigs
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/protein.sbt.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches:' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out)
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


@utils.in_tempdir
def test_sbt_hp_command_index(c):
    # test command-line creation of SBT database with hp sigs
    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/hp/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    db_out = c.output('hp.sbt.zip')

    c.run_sourmash('index', db_out, sigfile1, sigfile2,
                   '--scaled', '100', '-k', '19', '--hp')

    db2 = load_sbt_index(db_out)

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0, ignore_abundance=True,
                         do_containment=False, best_only=False)
    assert len(results) == 2

    results = db2.gather(sig2)
    assert results[0][0] == 1.0
    assert results[0][2] == db2._location
    assert results[0][2] == db_out


@utils.in_thisdir
def test_sbt_hp_command_search(c):
    # test command-line search/gather of SBT database with hp sigs
    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/hp.sbt.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches:' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


@utils.in_tempdir
def test_sbt_dayhoff_command_index(c):
    # test command-line creation of SBT database with dayhoff sigs
    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/dayhoff/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    db_out = c.output('dayhoff.sbt.zip')

    c.run_sourmash('index', db_out, sigfile1, sigfile2,
                   '--scaled', '100', '-k', '19', '--dayhoff')

    db2 = load_sbt_index(db_out)

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0, ignore_abundance=True,
                         do_containment=False, best_only=False)
    assert len(results) == 2

    results = db2.gather(sig2)
    assert results[0][0] == 1.0
    assert results[0][2] == db2._location
    assert results[0][2] == db_out


@utils.in_thisdir
def test_sbt_dayhoff_command_search(c):
    # test command-line search/gather of SBT database with dayhoff sigs
    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/dayhoff.sbt.zip')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches:' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


def test_sbt_node_cache():
    tree = SBT.load(utils.get_test_data('v6.sbt.json'),
                    leaf_loader=SigLeaf.load,
                    cache_size=1)

    testdata1 = utils.get_test_data(utils.SIG_FILES[0])
    to_search = load_one_signature(testdata1)

    # note: earlier versions of this test did containment on
    # the num MinHash in `to_search`, which doesn't work properly.
    # (See test_sbt_no_containment_on_num for test). So, to
    # fix the test for the new search API, we had to adjust
    # the threshold.
    search_obj = make_jaccard_search_query(threshold=0.05)
    results = list(tree.find(search_obj, to_search))
    assert len(results) == 4

    assert tree._nodescache.currsize == 1
    assert tree._nodescache.currsize == 1


def test_sbt_no_containment_on_num():
    tree = SBT.load(utils.get_test_data('v6.sbt.json'),
                    leaf_loader=SigLeaf.load,
                    cache_size=1)

    testdata1 = utils.get_test_data(utils.SIG_FILES[0])
    to_search = load_one_signature(testdata1)

    search_obj = make_jaccard_search_query(do_containment=True, threshold=0.05)
    with pytest.raises(TypeError) as exc:
        results = list(tree.find(search_obj, to_search))

    assert "this search requires a scaled signature" in str(exc)
