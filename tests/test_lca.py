"""
Tests for the 'sourmash lca' command line and high level API.
"""
import os
import shutil
import csv
import pytest
import glob

import sourmash_tst_utils as utils
import sourmash
from sourmash import load_one_signature, SourmashSignature

from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import LineagePair


def test_api_create_search():
    # create a database and then search for result.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    count = lca_db.insert(ss)
    assert count == len(ss.minhash)

    results = lca_db.search(ss, threshold=0.0)
    print(results)
    assert len(results) == 1
    (similarity, match, filename) = results[0]
    assert match.minhash == ss.minhash


def test_api_create_insert():
    # test some internal implementation stuff: create & then insert a sig.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss)

    ident = ss.name
    assert len(lca_db.ident_to_name) == 1
    assert ident in lca_db.ident_to_name
    assert lca_db.ident_to_name[ident] == ident
    assert len(lca_db.ident_to_idx) == 1
    assert lca_db.ident_to_idx[ident] == 0
    assert len(lca_db.hashval_to_idx) == len(ss.minhash)
    assert len(lca_db.idx_to_ident) == 1
    assert lca_db.idx_to_ident[0] == ident

    set_of_values = set()
    for vv in lca_db.hashval_to_idx.values():
        set_of_values.update(vv)
    assert len(set_of_values) == 1
    assert set_of_values == { 0 }

    assert not lca_db.idx_to_lid          # no lineage added
    assert not lca_db.lid_to_lineage      # no lineage added


def test_api_create_insert_bad_ksize():
    # can we insert a ksize=21 signature into a ksize=31 DB? hopefully not.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=21, scaled=1000)
    with pytest.raises(ValueError):
        lca_db.insert(ss)


def test_api_create_insert_bad_ident():
    # can we insert a signature with no/empty ident?
    ss1 = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                      ksize=31)
    ss2 = sourmash.load_one_signature(utils.get_test_data('63.fa.sig'),
                                      ksize=31)
    ss1.name = ''
    ss1.filename = ''
    ss2.name = ''
    ss2.filename = ''

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss1)
    lca_db.insert(ss2)
    # SUCCESS!
    # would fail, previously :)


def test_api_create_insert_bad_scaled():
    # can we insert a scaled=1000 signature into a scaled=500 DB?
    # hopefully not.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    assert ss.minhash.scaled == 1000

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=500)
    with pytest.raises(ValueError):
        lca_db.insert(ss)


def test_api_create_insert_bad_moltype():
    # can we insert a DNAsignature into a protein DB?
    # hopefully not.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    assert ss.minhash.moltype == 'DNA'

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=500, moltype='protein')
    with pytest.raises(ValueError):
        lca_db.insert(ss)


def test_api_create_insert_ident():
    # test some internal implementation stuff: signature inserted with
    # different ident than name.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss, ident='foo')

    ident = 'foo'
    assert len(lca_db.ident_to_name) == 1
    assert ident in lca_db.ident_to_name
    assert lca_db.ident_to_name[ident] == ss.name
    assert len(lca_db.ident_to_idx) == 1
    assert lca_db.ident_to_idx[ident] == 0
    assert len(lca_db.hashval_to_idx) == len(ss.minhash)
    assert len(lca_db.idx_to_ident) == 1
    assert lca_db.idx_to_ident[0] == ident

    set_of_values = set()
    for vv in lca_db.hashval_to_idx.values():
        set_of_values.update(vv)
    assert len(set_of_values) == 1
    assert set_of_values == { 0 }

    assert not lca_db.idx_to_lid          # no lineage added
    assert not lca_db.lid_to_lineage      # no lineage added
    assert not lca_db.lineage_to_lid
    assert not lca_db.lid_to_idx


def test_api_create_insert_two():
    # check internal details if multiple signatures are inserted.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    ss2 = sourmash.load_one_signature(utils.get_test_data('63.fa.sig'),
                                      ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss, ident='foo')
    lca_db.insert(ss2, ident='bar')

    ident = 'foo'
    ident2 = 'bar'
    assert len(lca_db.ident_to_name) == 2
    assert ident in lca_db.ident_to_name
    assert ident2 in lca_db.ident_to_name
    assert lca_db.ident_to_name[ident] == ss.name
    assert lca_db.ident_to_name[ident2] == ss2.name

    assert len(lca_db.ident_to_idx) == 2
    assert lca_db.ident_to_idx[ident] == 0
    assert lca_db.ident_to_idx[ident2] == 1

    combined_mins = set(ss.minhash.hashes.keys())
    combined_mins.update(set(ss2.minhash.hashes.keys()))
    assert len(lca_db.hashval_to_idx) == len(combined_mins)

    assert len(lca_db.idx_to_ident) == 2
    assert lca_db.idx_to_ident[0] == ident
    assert lca_db.idx_to_ident[1] == ident2

    set_of_values = set()
    for vv in lca_db.hashval_to_idx.values():
        set_of_values.update(vv)
    assert len(set_of_values) == 2
    assert set_of_values == { 0, 1 }

    assert not lca_db.idx_to_lid          # no lineage added
    assert not lca_db.lid_to_lineage      # no lineage added
    assert not lca_db.lineage_to_lid
    assert not lca_db.lid_to_idx


def test_api_create_insert_w_lineage():
    # test some internal implementation stuff - insert signature w/lineage
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lineage = ((LineagePair('rank1', 'name1'),
                LineagePair('rank2', 'name2')))

    lca_db.insert(ss, lineage=lineage)

    # basic ident stuff
    ident = ss.name
    assert len(lca_db.ident_to_name) == 1
    assert ident in lca_db.ident_to_name
    assert lca_db.ident_to_name[ident] == ident
    assert len(lca_db.ident_to_idx) == 1
    assert lca_db.ident_to_idx[ident] == 0
    assert len(lca_db.hashval_to_idx) == len(ss.minhash)
    assert len(lca_db.idx_to_ident) == 1
    assert lca_db.idx_to_ident[0] == ident

    # all hash values added
    set_of_values = set()
    for vv in lca_db.hashval_to_idx.values():
        set_of_values.update(vv)
    assert len(set_of_values) == 1
    assert set_of_values == { 0 }

    # check lineage stuff
    assert len(lca_db.idx_to_lid) == 1
    assert lca_db.idx_to_lid[0] == 0
    assert len(lca_db.lid_to_lineage) == 1
    assert lca_db.lid_to_lineage[0] == lineage
    assert lca_db.lid_to_idx[0] == { 0 }

    assert len(lca_db.lineage_to_lid) == 1
    assert lca_db.lineage_to_lid[lineage] == 0


def test_api_create_insert_w_bad_lineage():
    # test some internal implementation stuff - insert signature w/bad lineage
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lineage = ([LineagePair('rank1', 'name1'),
                LineagePair('rank2', 'name2')],)

    with pytest.raises(ValueError):
        lca_db.insert(ss, lineage=lineage)


def test_api_create_insert_w_bad_lineage_2():
    # test some internal implementation stuff - insert signature w/bad lineage
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lineage = 1 # something non-iterable...

    with pytest.raises(ValueError):
        lca_db.insert(ss, lineage=lineage)


def test_api_create_gather():
    # create a database, and then run gather on it.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss)

    results = lca_db.gather(ss, threshold_bp=0)
    print(results)
    assert len(results) == 1
    (similarity, match, filename) = results[0]
    assert match.minhash == ss.minhash


def test_api_add_genome_lineage():
    # LCA_Databases can store/retrieve arbitrary lineages/taxonomies.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    lineage = ((LineagePair('rank1', 'name1'),
               (LineagePair('rank2', 'name2'))))

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss, lineage=lineage)

    somehash = next(iter(ss.minhash.hashes.keys()))

    lineages = lca_db.get_lineage_assignments(somehash)
    assert len(lineages) == 1
    assert lineage in lineages


def test_api_insert_update():
    # check that cached parts of LCA_Database are updated when a new
    # signature is inserted.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    ss2 = sourmash.load_one_signature(utils.get_test_data('63.fa.sig'),
                                      ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss)

    all_mh = [ x.minhash for x in lca_db.signatures() ]
    assert ss.minhash in all_mh

    # see decorator @cached_property
    assert hasattr(lca_db, '_cache')
    assert lca_db._cache
    # inserting a signature should delete the cache
    lca_db.insert(ss2)
    assert not hasattr(lca_db, '_cache')

    # check that it's rebuilt etc. etc.
    all_mh = [ x.minhash for x in lca_db.signatures() ]
    assert ss.minhash in all_mh
    assert ss2.minhash in all_mh


def test_api_insert_retrieve_check_name():
    # check that signatures retrieved from LCA_Database objects have the
    # right name.
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss)

    sigs = list(lca_db.signatures())
    assert len(sigs) == 1
    retrieved_sig = sigs[0]
    assert retrieved_sig.name == ss.name
    assert retrieved_sig.minhash == ss.minhash


def test_api_create_insert_two_then_scale():
    # construct database, THEN downsample
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    ss2 = sourmash.load_one_signature(utils.get_test_data('63.fa.sig'),
                                      ksize=31)

    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    lca_db.insert(ss)
    lca_db.insert(ss2)

    # downsample everything to 5000
    lca_db.downsample_scaled(5000)

    ss.minhash = ss.minhash.downsample(scaled=5000)
    ss2.minhash = ss2.minhash.downsample(scaled=5000)

    # & check...
    combined_mins = set(ss.minhash.hashes.keys())
    combined_mins.update(set(ss2.minhash.hashes.keys()))
    assert len(lca_db.hashval_to_idx) == len(combined_mins)


def test_api_create_insert_scale_two():
    # downsample while constructing database
    ss = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'),
                                     ksize=31)
    ss2 = sourmash.load_one_signature(utils.get_test_data('63.fa.sig'),
                                      ksize=31)

    # downsample to 5000 while inserting:
    lca_db = sourmash.lca.LCA_Database(ksize=31, scaled=5000)
    count = lca_db.insert(ss)
    assert count == 1037
    assert count == len(ss.minhash.downsample(scaled=5000))
    lca_db.insert(ss2)

    # downsample sigs to 5000
    ss.minhash = ss.minhash.downsample(scaled=5000)
    ss2.minhash = ss2.minhash.downsample(scaled=5000)

    # & check...
    combined_mins = set(ss.minhash.hashes.keys())
    combined_mins.update(set(ss2.minhash.hashes.keys()))
    assert len(lca_db.hashval_to_idx) == len(combined_mins)


def test_load_single_db():
    filename = utils.get_test_data('lca/delmont-1.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    print(db)

    assert ksize == 31
    assert scaled == 10000


def test_databases():
    filename1 = utils.get_test_data('lca/delmont-1.lca.json')
    filename2 = utils.get_test_data('lca/delmont-2.lca.json')
    dblist, ksize, scaled = lca_utils.load_databases([filename1, filename2])

    print(dblist)

    assert len(dblist) == 2
    assert ksize == 31
    assert scaled == 10000


def test_databases_load_fail_on_no_JSON():
    filename1 = utils.get_test_data('prot/protein.zip')
    with pytest.raises(ValueError) as exc:
        dblist, ksize, scaled = lca_utils.load_databases([filename1])

    err = str(exc.value)
    print(err)
    assert f"'{filename1}' is not an LCA database file." in err


def test_databases_load_fail_on_dir():
    filename1 = utils.get_test_data('lca')
    with pytest.raises(ValueError) as exc:
        dblist, ksize, scaled = lca_utils.load_databases([filename1])

    err = str(exc.value)
    print(err)
    assert f"'{filename1}' is not a file and cannot be loaded as an LCA database" in err
    assert not 'found 0 matches total;' in err


def test_databases_load_fail_on_not_exist():
    filename1 = utils.get_test_data('does-not-exist')
    with pytest.raises(ValueError) as exc:
        dblist, ksize, scaled = lca_utils.load_databases([filename1])

    err = str(exc.value)
    print(err)
    assert f"'{filename1}' is not a file and cannot be loaded as an LCA database" in err
    assert not 'found 0 matches total;' in err

def test_db_repr():
    filename = utils.get_test_data('lca/delmont-1.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    assert repr(db) == "LCA_Database('{}')".format(filename)


def test_lca_index_signatures_method():
    # test 'signatures' method from base class Index
    filename = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    siglist = list(db.signatures())
    assert len(siglist) == 2


def test_lca_index_select():
    # test 'select' method from Index base class.

    filename = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(filename)

    xx = db.select(ksize=31)
    assert xx == db

    xx = db.select(moltype='DNA')
    assert xx == db

    with pytest.raises(ValueError):
        db.select(ksize=21)

    with pytest.raises(ValueError):
        db.select(moltype='protein')


def test_search_db_scaled_gt_sig_scaled():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)
    sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))

    results = db.search(sig, threshold=.01, ignore_abundance=True)
    match_sig = results[0][1]

    sig.minhash = sig.minhash.downsample(scaled=10000)
    assert sig.minhash == match_sig.minhash


def test_search_db_scaled_lt_sig_scaled():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)
    sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))
    sig.minhash = sig.minhash.downsample(scaled=100000)

    results = db.search(sig, threshold=.01, ignore_abundance=True)
    print(results)
    assert results[0][0] == 1.0
    match = results[0][1]

    orig_sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))
    assert orig_sig.minhash.jaccard(match.minhash, downsample=True) == 1.0


def test_gather_db_scaled_gt_sig_scaled():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)
    sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))

    results = db.gather(sig, threshold=.01, ignore_abundance=True)
    match_sig = results[0][1]

    sig.minhash = sig.minhash.downsample(scaled=10000)
    assert sig.minhash == match_sig.minhash


def test_gather_db_scaled_lt_sig_scaled():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)
    sig = sourmash.load_one_signature(utils.get_test_data('47.fa.sig'))
    sig.minhash = sig.minhash.downsample(scaled=100000)

    results = db.gather(sig, threshold=.01, ignore_abundance=True)
    match_sig = results[0][1]

    match_sig.minhash = match_sig.minhash.downsample(scaled=100000)
    assert sig.minhash == match_sig.minhash


def test_db_lineage_to_lid():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)

    d = db.lineage_to_lid
    items = list(d.items())
    items.sort()
    assert len(items) == 2

    print(items)

    lin1 = items[0][0][-1]
    assert lin1.rank == 'strain'
    assert lin1.name == 'Shewanella baltica OS185'
    lin1 = items[1][0][-1]
    assert lin1.rank == 'strain'
    assert lin1.name == 'Shewanella baltica OS223'


def test_db_lid_to_idx():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)

    d = db.lid_to_idx
    items = list(d.items())
    items.sort()
    assert len(items) == 2

    print(items)
    assert items == [(32, {32}), (48, {48})]


def test_db_idx_to_ident():
    dbfile = utils.get_test_data('lca/47+63.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(dbfile)

    d = db.idx_to_ident
    items = list(d.items())
    items.sort()
    assert len(items) == 2

    print(items)
    assert items == [(32, 'NC_009665'), (48, 'NC_011663')]


## command line tests


def test_run_sourmash_lca():
    status, out, err = utils.runscript('sourmash', ['lca'], fail_ok=True)
    assert status != 0                    # no args provided, ok ;)


def test_basic_index():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-1.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert 'Building LCA database with ksize=31 scaled=10000 moltype=DNA' in err
        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err


def test_basic_index_bad_spreadsheet():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/bad-spreadsheet.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err


def test_basic_index_broken_spreadsheet():
    # duplicate identifiers in this spreadsheet
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/bad-spreadsheet-2.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd, fail_ok=True)

        assert status != 0
        assert "multiple lineages for identifier TARA_ASE_MAG_00031" in err


def test_basic_index_too_many_strains_too_few_species():
    # explicit test for #841, where 'n_species' wasn't getting counted
    # if lineage was at strain level resolution.
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/podar-lineage.csv')
        input_sig = utils.get_test_data('47.fa.sig')
        lca_db = os.path.join(location, 'out.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig,
               '-C', '3', '--split-identifiers']
        status, out, err = utils.runscript('sourmash', cmd, fail_ok=True)

        assert not 'error: fewer than 20% of lineages' in err
        assert status == 0


def test_basic_index_too_few_species():
    # spreadsheets with too few species should be flagged, unless -f specified
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/tully-genome-sigs.classify.csv')

        # (these don't really matter, should break on load spreadsheet)
        input_sig = utils.get_test_data('47.fa.sig')
        lca_db = os.path.join(location, 'out.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '-C', '3']
        status, out, err = utils.runscript('sourmash', cmd, fail_ok=True)

        assert not '"ERROR: fewer than 20% of lineages have species-level resolution' in err
        assert status != 0


def test_basic_index_require_taxonomy():
    # no taxonomy in here
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/bad-spreadsheet-3.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', '--require-taxonomy', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd, fail_ok=True)

        assert status != 0
        assert "ERROR: no hash values found - are there any signatures?" in err


def test_basic_index_column_start():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-3.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', '-C', '3', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err


@utils.in_tempdir
def test_index_empty_sketch_name(c):
    # create two signatures with empty 'name' attributes
    cmd = ['sketch', 'dna', utils.get_test_data('genome-s12.fa.gz'),
           utils.get_test_data('genome-s11.fa.gz')]
    c.run_sourmash(*cmd)

    sig1 = c.output('genome-s11.fa.gz.sig')
    assert os.path.exists(sig1)
    sig2 = c.output('genome-s12.fa.gz.sig')
    assert os.path.exists(sig2)

    # can we insert them both?
    taxcsv = utils.get_test_data('lca/delmont-1.csv')
    cmd = ['lca', 'index', taxcsv, 'zzz', sig1, sig2]
    c.run_sourmash(*cmd)
    assert os.path.exists(c.output('zzz.lca.json'))

    print(c.last_result.out)
    print(c.last_result.err)
    assert 'WARNING: no lineage provided for 2 sig' in c.last_result.err


def test_basic_index_and_classify_with_tsv_and_gz():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-1.tsv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json.gz')

        cmd = ['lca', 'index', '--tabs', '--no-header', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_basic_index_and_classify():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-1.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_index_traverse():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-1.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        in_dir = os.path.join(location, 'sigs')
        os.mkdir(in_dir)
        shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

        cmd = ['lca', 'index', taxcsv, lca_db, in_dir]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err
        assert 'WARNING: 1 duplicate signatures.' not in err


@utils.in_tempdir
def test_index_traverse_force(c):
    # test the use of --force to load all files, not just .sig
    taxcsv = utils.get_test_data('lca/delmont-1.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = c.output('delmont-1.lca.json')

    in_dir = c.output('sigs')
    os.mkdir(in_dir)
    # name signature .txt instead of .sig:
    shutil.copyfile(input_sig, os.path.join(in_dir, 'q.txt'))

    # use --force
    cmd = ['lca', 'index', taxcsv, lca_db, in_dir, '-f']
    c.run_sourmash(*cmd)

    out = c.last_result.out
    err = c.last_result.err
    print(out)
    print(err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err
    assert 'WARNING: 1 duplicate signatures.' not in err


@utils.in_tempdir
def test_index_from_file_cmdline_sig(c):
    taxcsv = utils.get_test_data('lca/delmont-1.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = c.output('delmont-1.lca.json')

    file_list = c.output('sigs.list')
    with open(file_list, 'wt') as fp:
        print(input_sig, file=fp)

    cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '--from-file', file_list]
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err
    assert 'WARNING: 1 duplicate signatures.' in err


@utils.in_tempdir
def test_index_from_file(c):
    taxcsv = utils.get_test_data('lca/delmont-1.csv')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    lca_db = c.output('delmont-1.lca.json')

    file_list = c.output('sigs.list')
    with open(file_list, 'wt') as fp:
        print(input_sig, file=fp)

    cmd = ['lca', 'index', taxcsv, lca_db, '--from-file', file_list]
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
    assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err


@utils.in_tempdir
def test_index_fail_on_num(c):
    # lca index should yield a decent error message when attempted on 'num'
    sigfile = utils.get_test_data('num/63.fa.sig')
    taxcsv = utils.get_test_data('lca/podar-lineage.csv')

    with pytest.raises(ValueError):
        c.run_sourmash('lca', 'index', taxcsv, 'xxx.lca.json', sigfile, '-C', '3')

    err = c.last_result.err
    print(err)

    assert 'ERROR: cannot insert signature ' in err
    assert 'ERROR: cannot downsample signature; is it a scaled signature?' in err


def test_index_traverse_real_spreadsheet_no_report():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/tara-delmont-SuppTable3.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '-f']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 957 distinct identifiers in spreadsheet.' in err
        assert 'WARNING: no signatures for 956 spreadsheet rows.' in err
        assert 'WARNING: 105 unused lineages.' in err
        assert '(You can use --report to generate a detailed report.)' in err


def test_index_traverse_real_spreadsheet_report():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/tara-delmont-SuppTable3.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')
        report_loc = os.path.join(location, 'report.txt')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig, '--report',
               report_loc, '-f']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 957 distinct identifiers in spreadsheet.' in err
        assert 'WARNING: no signatures for 956 spreadsheet rows.' in err
        assert 'WARNING: 105 unused lineages.' in err
        assert '(You can use --report to generate a detailed report.)' not in err
        assert os.path.exists(report_loc)


def test_single_classify():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

        cmd = ['lca', 'classify', '--db', db1, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_single_classify_to_output():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

        cmd = ['lca', 'classify', '--db', db1, '--query', input_sig,
               '-o', os.path.join(location, 'outfile.txt')]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        outdata = open(os.path.join(location, 'outfile.txt'), 'rt').read()
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in outdata
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_single_classify_to_output_no_name():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        ss = sourmash.load_one_signature(input_sig, ksize=31)

        outsig_filename = os.path.join(location, 'q.sig')
        with open(outsig_filename, 'wt') as fp:
            # remove name from signature here --
            new_sig = sourmash.SourmashSignature(ss.minhash, filename='xyz')
            sourmash.save_signatures([new_sig], fp)

        cmd = ['lca', 'classify', '--db', db1, '--query', outsig_filename,
               '-o', os.path.join(location, 'outfile.txt')]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        outdata = open(os.path.join(location, 'outfile.txt'), 'rt').read()
        print((outdata,))
        assert 'xyz,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in outdata
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_single_classify_empty():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/both.lca.json')
        input_sig = utils.get_test_data('GCF_000005845.2_ASM584v2_genomic.fna.gz.sig')

        cmd = ['lca', 'classify', '--db', db1, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'GCF_000005845,nomatch,,,,,,,,' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_single_classify_traverse():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        in_dir = os.path.join(location, 'sigs')
        os.mkdir(in_dir)
        shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

        cmd = ['lca', 'classify', '--db', db1, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,Alteromonas,Alteromonas_macleodii' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_multi_query_classify_traverse():
    with utils.TempDirectory() as location:
        # both.lca.json is built from both dir and dir2
        db1 = utils.get_test_data('lca/both.lca.json')
        dir1 = utils.get_test_data('lca/dir1')
        dir2 = utils.get_test_data('lca/dir2')

        cmd = ['lca', 'classify', '--db', db1, '--query', dir1, dir2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        with open(utils.get_test_data('lca/classify-by-both.csv'), 'rt') as fp:
            fp_lines = fp.readlines()
            out_lines = out.splitlines()

            fp_lines.sort()
            out_lines.sort()

            assert len(fp_lines) == len(out_lines)
            for line1, line2 in zip(fp_lines, out_lines):
                assert line1.strip() == line2.strip(), (line1, line2)


@utils.in_tempdir
def test_multi_query_classify_query_from_file(c):
    # both.lca.json is built from both dir and dir2
    db1 = utils.get_test_data('lca/both.lca.json')
    dir1_glob = utils.get_test_data('lca/dir1/*.sig')
    dir1_files = glob.glob(dir1_glob)
    dir2_glob = utils.get_test_data('lca/dir2/*.sig')
    dir2_files = glob.glob(dir2_glob)

    file_list = c.output('file.list')
    with open(file_list, 'wt') as fp:
        print("\n".join(dir1_files), file=fp)
        print("\n".join(dir2_files), file=fp)

    cmd = ['lca', 'classify', '--db', db1, '--query-from-file', file_list]
    c.run_sourmash(*cmd)
    out = c.last_result.out

    with open(utils.get_test_data('lca/classify-by-both.csv'), 'rt') as fp:
        fp_lines = fp.readlines()
        out_lines = out.splitlines()

        fp_lines.sort()
        out_lines.sort()

        assert len(fp_lines) == len(out_lines)
        for line1, line2 in zip(fp_lines, out_lines):
            assert line1.strip() == line2.strip(), (line1, line2)


@utils.in_tempdir
def test_multi_query_classify_query_from_file_and_query(c):
    # both.lca.json is built from both dir and dir2
    db1 = utils.get_test_data('lca/both.lca.json')
    dir1_glob = utils.get_test_data('lca/dir1/*.sig')
    dir1_files = glob.glob(dir1_glob)
    dir2_glob = utils.get_test_data('lca/dir2/*.sig')
    dir2_files = glob.glob(dir2_glob)

    file_list = c.output('file.list')
    with open(file_list, 'wt') as fp:
        print("\n".join(dir1_files[1:]), file=fp)   # leave off first one
        print("\n".join(dir2_files), file=fp)

    cmd = ['lca', 'classify', '--db', db1, '--query', dir1_files[0],
           '--query-from-file', file_list]
    c.run_sourmash(*cmd)
    out = c.last_result.out

    with open(utils.get_test_data('lca/classify-by-both.csv'), 'rt') as fp:
        fp_lines = fp.readlines()
        out_lines = out.splitlines()

        fp_lines.sort()
        out_lines.sort()

        assert len(fp_lines) == len(out_lines)
        for line1, line2 in zip(fp_lines, out_lines):
            assert line1.strip() == line2.strip(), (line1, line2)


def test_multi_db_multi_query_classify_traverse():
    with utils.TempDirectory() as location:
        # two halves of both.lca.json, see above test.
        db1 = utils.get_test_data('lca/dir1.lca.json')
        db2 = utils.get_test_data('lca/dir2.lca.json')
        dir1 = utils.get_test_data('lca/dir1')
        dir2 = utils.get_test_data('lca/dir2')

        cmd = ['lca', 'classify', '--db', db1, db2, '--query', dir1, dir2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        with open(utils.get_test_data('lca/classify-by-both.csv'), 'rt') as fp:
            fp_lines = fp.readlines()
            out_lines = out.splitlines()

            fp_lines.sort()
            out_lines.sort()

            assert len(fp_lines) == len(out_lines)
            for line1, line2 in zip(fp_lines, out_lines):
                assert line1.strip() == line2.strip(), (line1, line2)


def test_unassigned_internal_index_and_classify():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-4.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,unassigned,Alteromonadaceae,unassigned,Alteromonas_macleodii' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_unassigned_last_index_and_classify():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-5.csv')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '1 identifiers used out of 1 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Alteromonadaceae,,,\r\n' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


def test_index_and_classify_internal_unassigned_multi():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-6.csv')
        input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig1, input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in err

        # classify input_sig1
        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig1]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,unassigned,unassigned,Alteromonadaceae,,,\r\n' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err

        # classify input_sig2
        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_PSW_MAG_00136,found,Eukaryota,Chlorophyta,Prasinophyceae,unassigned,unassigned,Ostreococcus,,\r\n' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 1 LCA databases' in err


@utils.in_tempdir
def test_classify_majority_vote_1(c):
    # classify merged signature using lca should yield no results

    # build database
    taxcsv = utils.get_test_data('lca/delmont-6.csv')
    input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
    lca_db = c.output('delmont-1.lca.json')

    c.run_sourmash('lca', 'index', taxcsv, lca_db, input_sig1, input_sig2)

    print(c.last_command)
    print(c.last_result.out)
    print(c.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in c.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in c.last_result.err
    assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in c.last_result.err

    # merge input_sig1 and input_sig2
    c.run_sourmash('signature', 'merge', input_sig1, input_sig2, '-k', '31', '--flatten', '-o', 'sig1and2.sig')
    sig1and2 = c.output('sig1and2.sig')

    # lca classify should yield no results
    c.run_sourmash('lca', 'classify', '--db', lca_db, '--query', sig1and2)

    print(c.last_command)
    print(c.last_result.out)
    print(c.last_result.err)

    assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in c.last_result.out
    assert 'disagree,,,,,,,,' in c.last_result.out
    assert 'classified 1 signatures total' in c.last_result.err
    assert 'loaded 1 LCA databases' in c.last_result.err



@utils.in_tempdir
def test_classify_majority_vote_2(c):
    # classify same signature with same database using --majority
    # should yield results

    # build database
    taxcsv = utils.get_test_data('lca/delmont-6.csv')
    input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
    lca_db = c.output('delmont-1.lca.json')

    c.run_sourmash('lca', 'index', taxcsv, lca_db, input_sig1, input_sig2)

    print(c.last_command)
    print(c.last_result.out)
    print(c.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in c.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in c.last_result.err
    assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in c.last_result.err

    # merge input_sig1 and input_sig2
    c.run_sourmash('signature', 'merge', input_sig1, input_sig2, '-k', '31', '--flatten', '-o', 'sig1and2.sig')
    sig1and2 = c.output('sig1and2.sig')

    # majority vote classify
    c.run_sourmash('lca', 'classify', '--db', lca_db, '--query', sig1and2, '--majority')

    print(c.last_command)
    print(c.last_result.out)
    print(c.last_result.err)

    assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in c.last_result.out
    assert 'found,Eukaryota,Chlorophyta,Prasinophyceae,unassigned,unassigned,Ostreococcus' in c.last_result.out
    assert 'classified 1 signatures total' in c.last_result.err
    assert 'loaded 1 LCA databases' in c.last_result.err


@utils.in_tempdir
def test_classify_majority_vote_3(c):
    # classify signature with nothing in counts

    # build database
    taxcsv = utils.get_test_data('lca/delmont-6.csv')
    input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
    lca_db = c.output('delmont-1.lca.json')

    c.run_sourmash('lca', 'index', taxcsv, lca_db, input_sig1, input_sig2)

    print(c.last_command)
    print(c.last_result.out)
    print(c.last_result.err)

    assert os.path.exists(lca_db)

    assert "** assuming column 'MAGs' is identifiers in spreadsheet" in c.last_result.err
    assert "** assuming column 'Domain' is superkingdom in spreadsheet" in c.last_result.err
    assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in c.last_result.err

    # obtain testdata '47.fa.sig'
    testdata1 = utils.get_test_data('47.fa.sig')

    # majority vote classify
    c.run_sourmash('lca', 'classify', '--db', lca_db, '--query', testdata1, '--majority')

    print(c.last_command)
    print(c.last_result.out)
    print(c.last_result.err)

    assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in c.last_result.out
    assert 'nomatch,,,,,,,,' in c.last_result.out
    assert 'classified 1 signatures total' in c.last_result.err
    assert 'loaded 1 LCA databases' in c.last_result.err


def test_multi_db_classify():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        db2 = utils.get_test_data('lca/delmont-2.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

        cmd = ['lca', 'classify', '--db', db1, db2, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'ID,status,superkingdom,phylum,class,order,family,genus,species' in out
        assert 'TARA_ASE_MAG_00031,found,Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,,,,' in out
        assert 'classified 1 signatures total' in err
        assert 'loaded 2 LCA databases' in err


def test_classify_unknown_hashes():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca-root/tax.csv')
        input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
        input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
        lca_db = os.path.join(location, 'lca-root.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert '1 identifiers used out of 2 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'classify', '--db', lca_db, '--query', input_sig1]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert '(root)' not in out
        assert 'TARA_MED_MAG_00029,found,Archaea,Euryarcheoata,unassigned,unassigned,novelFamily_I' in out


def test_single_summarize():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

        cmd = ['lca', 'summarize', '--db', db1, '--query', input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'loaded 1 signatures from 1 files total.' in err
        assert '100.0%   200   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales' in out


def test_single_summarize_singleton():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')

        cmd = ['lca', 'summarize', '--db', db1, '--query', input_sig,]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'loaded 1 signatures from 1 files total.' in err
        assert '100.0%   200   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales' in out
        assert 'test-data/lca/TARA_ASE_MAG_00031.sig:5b438c6c TARA_ASE_MAG_00031' in out


@utils.in_tempdir
def test_single_summarize_traverse(c):
    db1 = utils.get_test_data('lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    in_dir = c.output('sigs')
    os.mkdir(in_dir)
    shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

    cmd = ['lca', 'summarize', '--db', db1, '--query', in_dir]
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert 'loaded 1 signatures from 1 files total.' in err
    assert '100.0%   200   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales' in out

@utils.in_tempdir
def test_single_summarize_singleton_traverse(c):
    db1 = utils.get_test_data('lca/delmont-1.lca.json')
    input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
    in_dir = c.output('sigs')
    os.mkdir(in_dir)
    shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

    cmd = ['lca', 'summarize', '--db', db1, '--query', in_dir]
    c.run_sourmash(*cmd)

    out = c.last_result.out
    print(out)
    err = c.last_result.err
    print(err)

    assert 'loaded 1 signatures from 1 files total.' in err
    assert '100.0%   200   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales' in out
    assert 'q.sig:5b438c6c TARA_ASE_MAG_00031' in out


def test_single_summarize_to_output():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        in_dir = os.path.join(location, 'sigs')
        os.mkdir(in_dir)
        shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

        cmd = ['lca', 'summarize', '--db', db1, '--query', input_sig,
               '-o', os.path.join(location, 'output.txt')]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        outdata = open(os.path.join(location, 'output.txt'), 'rt').read()

        assert 'loaded 1 signatures from 1 files total.' in err
        assert '200,Bacteria,Proteobacteria,Gammaproteobacteria' in outdata


def test_single_summarize_scaled():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/delmont-1.lca.json')
        input_sig = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        in_dir = os.path.join(location, 'sigs')
        os.mkdir(in_dir)
        shutil.copyfile(input_sig, os.path.join(in_dir, 'q.sig'))

        cmd = ['lca', 'summarize', '--db', db1, '--query', input_sig,
               '--scaled', '100000']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'loaded 1 signatures from 1 files total.' in err
        assert '100.0%    27   Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales'


def test_multi_summarize_with_unassigned_singleton():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca/delmont-6.csv')
        input_sig1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.sig')
        input_sig2 = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig1, input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig1,
               input_sig2, '--ignore-abundance']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'loaded 2 signatures from 2 files total.' in err

        out_lines = out.splitlines()
        def remove_line_startswith(x, check=None):
           for line in out_lines:
               if line.startswith(x):
                   out_lines.remove(line)
                   if check:
                       # make sure the check value is in there
                       assert check in line
                   return line
           assert 0, "couldn't find {}".format(x)

        # note, proportions/percentages are now per-file
        remove_line_startswith('100.0%   200   Bacteria ', 'TARA_ASE_MAG_00031.sig:5b438c6c')
        remove_line_startswith('100.0%   200   Bacteria;Proteobacteria;unassigned;unassigned ')
        remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta ')
        remove_line_startswith('100.0%  1231   Eukaryota ', 'TARA_PSW_MAG_00136.sig:db50b713')
        remove_line_startswith('100.0%   200   Bacteria;Proteobacteria ')
        remove_line_startswith('100.0%   200   Bacteria;Proteobacteria;unassigned ')
        remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae ')
        remove_line_startswith('100.0%   200   Bacteria;Proteobacteria;unassigned;unassigned;Alteromonadaceae ')
        remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae;unassigned;unassigned ')
        remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae;unassigned ')
        remove_line_startswith('100.0%  1231   Eukaryota;Chlorophyta;Prasinophyceae;unassigned;unassigned;Ostreococcus ')
        assert not out_lines


def test_summarize_to_root():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca-root/tax.csv')
        input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
        input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
        lca_db = os.path.join(location, 'lca-root.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig1, input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig2,
               '--ignore-abundance']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert '78.6%    99   Archaea' in out
        assert '21.4%    27   (root)' in out


def test_summarize_unknown_hashes():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca-root/tax.csv')
        input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
        input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
        lca_db = os.path.join(location, 'lca-root.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert '1 identifiers used out of 2 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig1]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert '(root)' not in out
        assert '11.5%    27   Archaea;Euryarcheoata;unassigned;unassigned;novelFamily_I' in out


def test_summarize_to_root_abund():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca-root/tax.csv')
        input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
        input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
        lca_db = os.path.join(location, 'lca-root.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig1, input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert '2 identifiers used out of 2 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert '78.9%   101   Archaea' in out
        assert '21.1%    27   (root)' in out


def test_summarize_unknown_hashes_abund():
    with utils.TempDirectory() as location:
        taxcsv = utils.get_test_data('lca-root/tax.csv')
        input_sig1 = utils.get_test_data('lca-root/TARA_MED_MAG_00029.fa.sig')
        input_sig2 = utils.get_test_data('lca-root/TOBG_MED-875.fna.gz.sig')
        lca_db = os.path.join(location, 'lca-root.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert '1 identifiers used out of 2 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'summarize', '--db', lca_db, '--query', input_sig1]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert '(root)' not in out
        assert '11.5%    27   Archaea;Euryarcheoata;unassigned;unassigned;novelFamily_I' in out


@utils.in_thisdir
def test_lca_summarize_abund_hmp(c):
    # test lca summarize --with-abundance on some real data
    queryfile = utils.get_test_data('hmp-sigs/G36354.sig.gz')
    dbname = utils.get_test_data('hmp-sigs/G36354-matches.lca.json.gz')

    c.run_sourmash('lca', 'summarize', '--db', dbname, '--query', queryfile)

    assert '32.1%  1080   p__Firmicutes;c__Bacilli;o__Lactobacillales' in c.last_result.out


@utils.in_thisdir
def test_lca_summarize_abund_fake_no_abund(c):
    # test lca summarize on some known/fake data; see docs for explanation.
    queryfile = utils.get_test_data('fake-abund/query.sig.gz')
    dbname = utils.get_test_data('fake-abund/matches.lca.json.gz')

    c.run_sourmash('lca', 'summarize', '--db', dbname, '--query', queryfile,
                   '--ignore-abundance')

    assert 'NOTE: discarding abundances in query, since --ignore-abundance' in c.last_result.err
    assert '79.6%   550   Bacteria' in c.last_result.out
    assert '20.4%   141   Archaea' in c.last_result.out


@utils.in_thisdir
def test_lca_summarize_abund_fake_yes_abund(c):
    # test lca summarize abundance weighting on some known/fake data
    queryfile = utils.get_test_data('fake-abund/query.sig.gz')
    dbname = utils.get_test_data('fake-abund/matches.lca.json.gz')

    c.run_sourmash('lca', 'summarize', '--db', dbname, '--query', queryfile)

    assert '43.2%   563   Bacteria' in c.last_result.out
    assert '56.8%   740   Archaea' in c.last_result.out


def test_rankinfo_on_multi():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/dir1.lca.json')
        db2 = utils.get_test_data('lca/dir2.lca.json')

        cmd = ['lca', 'rankinfo', db1, db2]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        lines = out.splitlines()
        lines.remove('superkingdom: 0 (0.0%)')
        lines.remove('phylum: 464 (12.8%)')
        lines.remove('class: 533 (14.7%)')
        lines.remove('order: 1050 (29.0%)')
        lines.remove('family: 695 (19.2%)')
        lines.remove('genus: 681 (18.8%)')
        lines.remove('species: 200 (5.5%)')
        lines.remove('strain: 0 (0.0%)')

        assert not lines


def test_rankinfo_on_single():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/both.lca.json')

        cmd = ['lca', 'rankinfo', db1]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        lines = out.splitlines()
        lines.remove('superkingdom: 0 (0.0%)')
        lines.remove('phylum: 464 (12.8%)')
        lines.remove('class: 533 (14.7%)')
        lines.remove('order: 1050 (29.0%)')
        lines.remove('family: 695 (19.2%)')
        lines.remove('genus: 681 (18.8%)')
        lines.remove('species: 200 (5.5%)')
        lines.remove('strain: 0 (0.0%)')

        assert not lines


def test_rankinfo_no_tax():
    with utils.TempDirectory() as location:
        # note: TARA_PSW_MAG_00136 is _not_ in delmont-1.csv.
        taxcsv = utils.get_test_data('lca/delmont-1.csv')
        input_sig = utils.get_test_data('lca/TARA_PSW_MAG_00136.sig')
        lca_db = os.path.join(location, 'delmont-1.lca.json')

        cmd = ['lca', 'index', taxcsv, lca_db, input_sig]
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert os.path.exists(lca_db)

        assert "** assuming column 'MAGs' is identifiers in spreadsheet" in err
        assert "** assuming column 'Domain' is superkingdom in spreadsheet" in err
        assert '0 identifiers used out of 1 distinct identifiers in spreadsheet.' in err

        cmd = ['lca', 'rankinfo', lca_db]
        status, out, err = utils.runscript('sourmash', cmd)


def test_rankinfo_with_min():
    with utils.TempDirectory() as location:
        db1 = utils.get_test_data('lca/dir1.lca.json')
        db2 = utils.get_test_data('lca/dir2.lca.json')

        cmd = ['lca', 'rankinfo', db1, db2, '--minimum-num', '1']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        lines = out.splitlines()
        lines.remove('superkingdom: 0 (0.0%)')
        lines.remove('phylum: 464 (12.8%)')
        lines.remove('class: 533 (14.7%)')
        lines.remove('order: 1050 (29.0%)')
        lines.remove('family: 695 (19.2%)')
        lines.remove('genus: 681 (18.8%)')
        lines.remove('species: 200 (5.5%)')
        lines.remove('strain: 0 (0.0%)')

        assert not lines


def test_compare_csv():
    with utils.TempDirectory() as location:
        a = utils.get_test_data('lca/classify-by-both.csv')
        b = utils.get_test_data('lca/tara-delmont-SuppTable3.csv')

        cmd = ['lca', 'compare_csv', a, b, '-f']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'loaded 106 distinct lineages, 957 rows' in err
        assert 'missing 937 assignments in classify spreadsheet.' in err
        assert '20 total assignments, 0 differ between spreadsheets.' in err


def test_compare_csv_real():
    with utils.TempDirectory() as location:
        a = utils.get_test_data('lca/tully-genome-sigs.classify.csv')
        b = utils.get_test_data('lca/tully-query.delmont-db.sigs.classify.csv')

        cmd = ['lca', 'compare_csv', a, b, '--start-column=3', '-f']
        status, out, err = utils.runscript('sourmash', cmd)

        print(cmd)
        print(out)
        print(err)

        assert 'loaded 87 distinct lineages, 2631 rows' in err
        assert 'missing 71 assignments in classify spreadsheet.' in err
        assert 'missing 1380 assignments in custom spreadsheet.' in err
        assert '(these will not be evaluated any further)' in err
        assert '987 total assignments, 889 differ between spreadsheets.' in err
        assert '296 are compatible (one lineage is ancestor of another.' in err
        assert '593 are incompatible (there is a disagreement in the trees).' in err
        assert '164 incompatible at rank superkingdom' in err
        assert '255 incompatible at rank phylum' in err
        assert '107 incompatible at rank class' in err
        assert '54 incompatible at rank order' in err
        assert '13 incompatible at rank family' in err
        assert '0 incompatible at rank genus' in err
        assert '0 incompatible at rank species' in err


@utils.in_tempdir
def test_incompat_lca_db_ksize_2(c):
    # test on gather - create a database with ksize of 25
    testdata1 = utils.get_test_data('lca/TARA_ASE_MAG_00031.fa.gz')
    c.run_sourmash('compute', '-k', '25', '--scaled', '1000', testdata1,
                   '-o', 'test_db.sig')
    print(c)

    c.run_sourmash('lca', 'index', utils.get_test_data('lca/delmont-1.csv',),
                   'test.lca.json', 'test_db.sig',
                    '-k', '25', '--scaled', '10000')
    print(c)

    # this should fail: the LCA database has ksize 25, and the query sig has
    # no compatible ksizes.
    with pytest.raises(ValueError) as e:
        c.run_sourmash('gather', utils.get_test_data('lca/TARA_ASE_MAG_00031.sig'), 'test.lca.json')

    err = c.last_result.err
    print(err)

    assert "ERROR: cannot use 'test.lca.json' for this query." in err
    assert "ksize on this database is 25; this is different from requested ksize of 31"


@utils.in_tempdir
def test_lca_index_empty(c):
    # test lca index with an empty taxonomy CSV, followed by a load & gather.
    sig2file = utils.get_test_data('2.fa.sig')
    sig47file = utils.get_test_data('47.fa.sig')
    sig63file = utils.get_test_data('63.fa.sig')

    sig63 = load_one_signature(sig63file, ksize=31)

    # create an empty spreadsheet
    with open(c.output('empty.csv'), 'wt') as fp:
        fp.write('accession,superkingdom,phylum,class,order,family,genus,species,strain')

    # index!
    c.run_sourmash('lca', 'index', 'empty.csv', 'xxx.lca.json',
                   sig2file, sig47file, sig63file, '--scaled', '1000')

    # can we load and search?
    lca_db_filename = c.output('xxx.lca.json')
    db, ksize, scaled = lca_utils.load_single_database(lca_db_filename)

    results = db.gather(sig63)
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig.minhash == sig63.minhash
    assert name == lca_db_filename


def test_lca_gather_threshold_1():
    # test gather() method, in some detail; see same tests for sbt.
    sig2file = utils.get_test_data('2.fa.sig')
    sig47file = utils.get_test_data('47.fa.sig')
    sig63file = utils.get_test_data('63.fa.sig')

    sig2 = load_one_signature(sig2file, ksize=31)
    sig47 = load_one_signature(sig47file, ksize=31)
    sig63 = load_one_signature(sig63file, ksize=31)

    # construct LCA Database
    db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    db.insert(sig2)
    db.insert(sig47)
    db.insert(sig63)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures all have scaled=1000.

    mins = list(sorted(sig2.minhash.hashes.keys()))
    new_mh = sig2.minhash.copy_and_clear()

    # query with empty hashes
    assert not new_mh
    with pytest.raises(ValueError):
        db.gather(SourmashSignature(new_mh))

    # add one hash
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 1

    results = db.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig.minhash == sig2.minhash
    assert name == None

    # check with a threshold -> should be no results.
    with pytest.raises(ValueError):
        db.gather(SourmashSignature(new_mh), threshold_bp=5000)

    # add three more hashes => length of 4
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    new_mh.add_hash(mins.pop())
    assert len(new_mh) == 4

    results = db.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig.minhash == sig2.minhash
    assert name == None

    # check with a too-high threshold -> should be no results.
    with pytest.raises(ValueError):
        db.gather(SourmashSignature(new_mh), threshold_bp=5000)


def test_lca_gather_threshold_5():
    # test gather() method, in some detail; see same tests for sbt.
    sig2file = utils.get_test_data('2.fa.sig')
    sig47file = utils.get_test_data('47.fa.sig')
    sig63file = utils.get_test_data('63.fa.sig')

    sig2 = load_one_signature(sig2file, ksize=31)
    sig47 = load_one_signature(sig47file, ksize=31)
    sig63 = load_one_signature(sig63file, ksize=31)

    # construct LCA Database
    db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    db.insert(sig2)
    db.insert(sig47)
    db.insert(sig63)

    # now construct query signatures with specific numbers of hashes --
    # note, these signatures both have scaled=1000.

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
    results = db.gather(SourmashSignature(new_mh))
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig.minhash == sig2.minhash
    assert name == None

    # now, check with a threshold_bp that should be meet-able.
    results = db.gather(SourmashSignature(new_mh), threshold_bp=5000)
    assert len(results) == 1
    containment, match_sig, name = results[0]
    assert containment == 1.0
    assert match_sig.minhash == sig2.minhash
    assert name == None


def test_gather_multiple_return():
    sig2file = utils.get_test_data('2.fa.sig')
    sig47file = utils.get_test_data('47.fa.sig')
    sig63file = utils.get_test_data('63.fa.sig')

    sig2 = load_one_signature(sig2file, ksize=31)
    sig47 = load_one_signature(sig47file, ksize=31)
    sig63 = load_one_signature(sig63file, ksize=31)

    # construct LCA Database
    db = sourmash.lca.LCA_Database(ksize=31, scaled=1000)
    db.insert(sig2)
    db.insert(sig47)
    db.insert(sig63)

    # now, run gather. how many results do we get, and are they in the
    # right order?
    results = db.gather(sig63)
    print(len(results))
    assert len(results) == 1
    assert results[0][0] == 1.0


def test_lca_db_protein_build():
    # test programmatic creation of LCA database with protein sigs in it
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    db = sourmash.lca.LCA_Database(ksize=19, scaled=100, moltype='protein')
    assert db.insert(sig1)
    assert db.insert(sig2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db.search(sig1, threshold=0.0)
    assert len(results) == 2

    results = db.gather(sig2)
    assert results[0][0] == 1.0


@utils.in_tempdir
def test_lca_db_protein_save_load(c):
    # test save/load of programmatically created db with protein sigs in it
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    db = sourmash.lca.LCA_Database(ksize=19, scaled=100, moltype='protein')
    assert db.insert(sig1)
    assert db.insert(sig2)

    db.save(c.output('xxx.lca.json'))
    del db

    x = sourmash.lca.lca_db.load_single_database(c.output('xxx.lca.json'))
    db2 = x[0]
    assert db2.moltype == 'protein'

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    print('XXX', mh_list[0].ksize)
    print('YYY', sig1.minhash.ksize)
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0)
    assert len(results) == 2

    results = db2.gather(sig2)
    assert results[0][0] == 1.0


@utils.in_tempdir
def test_lca_db_protein_command_index(c):
    # test command-line creation of LCA database with protein sigs
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')
    lineages = utils.get_test_data('prot/gtdb-subset-lineages.csv')

    db_out = c.output('protein.lca.json')

    c.run_sourmash('lca', 'index', lineages, db_out, sigfile1, sigfile2,
                   '-C', '3', '--split-identifiers', '--require-taxonomy',
                   '--scaled', '100', '-k', '19', '--protein')

    x = sourmash.lca.lca_db.load_single_database(db_out)
    db2 = x[0]
    assert db2.moltype == 'protein'

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0)
    assert len(results) == 2

    results = db2.gather(sig2)
    assert results[0][0] == 1.0


@utils.in_thisdir
def test_lca_db_protein_command_search(c):
    # test command-line search/gather of LCA database with protein sigs
    # (LCA database created as above)
    sigfile1 = utils.get_test_data('prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/protein.lca.json.gz')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches:' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out)
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


def test_lca_db_hp_build():
    # test programmatic creation of LCA database with hp sigs in it
    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/hp/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    db = sourmash.lca.LCA_Database(ksize=19, scaled=100, moltype='hp')
    assert db.insert(sig1)
    assert db.insert(sig2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db.search(sig1, threshold=0.0)
    assert len(results) == 2

    results = db.gather(sig2)
    assert results[0][0] == 1.0


@utils.in_tempdir
def test_lca_db_hp_save_load(c):
    # test save/load of programmatically created db with hp sigs in it
    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/hp/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    db = sourmash.lca.LCA_Database(ksize=19, scaled=100, moltype='hp')
    assert db.insert(sig1)
    assert db.insert(sig2)

    db.save(c.output('xxx.lca.json'))
    del db

    x = sourmash.lca.lca_db.load_single_database(c.output('xxx.lca.json'))
    db2 = x[0]
    assert db2.moltype == 'hp'

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0)
    assert len(results) == 2

    results = db2.gather(sig2)
    assert results[0][0] == 1.0


@utils.in_tempdir
def test_lca_db_hp_command_index(c):
    # test command-line creation of LCA database with hp sigs
    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/hp/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')
    lineages = utils.get_test_data('prot/gtdb-subset-lineages.csv')

    db_out = c.output('hp.lca.json')

    c.run_sourmash('lca', 'index', lineages, db_out, sigfile1, sigfile2,
                   '-C', '3', '--split-identifiers', '--require-taxonomy',
                   '--scaled', '100', '-k', '19', '--hp')

    x = sourmash.lca.lca_db.load_single_database(db_out)
    db2 = x[0]
    assert db2.moltype == 'hp'

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0)
    assert len(results) == 2

    results = db2.gather(sig2)
    assert results[0][0] == 1.0


@utils.in_thisdir
def test_lca_db_hp_command_search(c):
    # test command-line search/gather of LCA database with hp sigs
    # (LCA database created as above)
    sigfile1 = utils.get_test_data('prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/hp.lca.json.gz')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches:' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out


def test_lca_db_dayhoff_build():
    # test programmatic creation of LCA database with dayhoff sigs in it
    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/dayhoff/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    db = sourmash.lca.LCA_Database(ksize=19, scaled=100, moltype='dayhoff')
    assert db.insert(sig1)
    assert db.insert(sig2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db.search(sig1, threshold=0.0)
    assert len(results) == 2

    results = db.gather(sig2)
    assert results[0][0] == 1.0


@utils.in_tempdir
def test_lca_db_dayhoff_save_load(c):
    # test save/load of programmatically created db with dayhoff sigs in it
    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/dayhoff/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    db = sourmash.lca.LCA_Database(ksize=19, scaled=100, moltype='dayhoff')
    assert db.insert(sig1)
    assert db.insert(sig2)

    db.save(c.output('xxx.lca.json'))
    del db

    x = sourmash.lca.lca_db.load_single_database(c.output('xxx.lca.json'))
    db2 = x[0]
    assert db2.moltype == 'dayhoff'

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0)
    assert len(results) == 2

    results = db2.gather(sig2)
    assert results[0][0] == 1.0


@utils.in_tempdir
def test_lca_db_dayhoff_command_index(c):
    # test command-line creation of LCA database with dayhoff sigs
    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    sigfile2 = utils.get_test_data('prot/dayhoff/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig')
    lineages = utils.get_test_data('prot/gtdb-subset-lineages.csv')

    db_out = c.output('dayhoff.lca.json')

    c.run_sourmash('lca', 'index', lineages, db_out, sigfile1, sigfile2,
                   '-C', '3', '--split-identifiers', '--require-taxonomy',
                   '--scaled', '100', '-k', '19', '--dayhoff')

    x = sourmash.lca.lca_db.load_single_database(db_out)
    db2 = x[0]
    assert db2.moltype == 'dayhoff'

    sig1 = sourmash.load_one_signature(sigfile1)
    sig2 = sourmash.load_one_signature(sigfile2)

    # check reconstruction --
    mh_list = [ x.minhash for x in db2.signatures() ]
    assert len(mh_list) == 2
    assert sig1.minhash in mh_list
    assert sig2.minhash in mh_list

    # and search, gather
    results = db2.search(sig1, threshold=0.0)
    assert len(results) == 2

    results = db2.gather(sig2)
    assert results[0][0] == 1.0


@utils.in_thisdir
def test_lca_db_dayhoff_command_search(c):
    # test command-line search/gather of LCA database with dayhoff sigs
    # (LCA database created as above)
    sigfile1 = utils.get_test_data('prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig')
    db_out = utils.get_test_data('prot/dayhoff.lca.json.gz')

    c.run_sourmash('search', sigfile1, db_out, '--threshold', '0.0')
    assert '2 matches:' in c.last_result.out

    c.run_sourmash('gather', sigfile1, db_out, '--threshold', '0.0')
    assert 'found 1 matches total' in c.last_result.out
    assert 'the recovered matches hit 100.0% of the query' in c.last_result.out
