from __future__ import print_function, unicode_literals

import os
import sourmash
from sourmash.index import LinearIndex
from sourmash_lib.sbt import SBT, GraphFactory, Leaf
from . import sourmash_tst_utils as utils


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
    print([s[1].name() for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2

    # search for sig47 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss47, threshold=0.1)
    print([s[1].name() for s in sr])
    assert len(sr) == 2
    sr.sort(key=lambda x: -x[0])
    assert sr[0][1] == ss47
    assert sr[1][1] == ss63

    # search for sig63 with lower threshold; search order not guaranteed.
    sr = lidx.search(ss63, threshold=0.1)
    print([s[1].name() for s in sr])
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
    assert len(matches) == 2
    assert matches[0][0] == 1.0
    assert matches[0][1] == ss47
    assert round(matches[1][0], 2) == 0.49
    assert matches[1][1] == ss63


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

        from sourmash import load_signatures
        si = set(load_signatures(filename))

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
    print([s[1].name() for s in sr])
    assert len(sr) == 1
    assert sr[0][1] == ss2


def test_linear_index_multik_select():
    # this loads three ksizes, 21/31/51
    sig2 = utils.get_test_data('2.fa.sig')
    siglist = sourmash.load_signatures(sig2)

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
    # this loads two ksizes(21, 30), and two moltypes (DNA and protein)
    filename = utils.get_test_data('genome-s10+s11.sig')
    siglist = sourmash.load_signatures(filename)

    linear = LinearIndex()
    for ss in siglist:
        linear.insert(ss)

    # select most specific DNA
    linear2 = linear.select(ksize=30, moltype='DNA')
    assert len(linear2) == 1

    # select most specific protein
    linear2 = linear.select(ksize=30, moltype='protein')
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
