"Tests for search.py code."

# CTB TODO: test search protocol with mock class?

import pytest

from sourmash import search, SourmashSignature, MinHash
from sourmash.search import make_jaccard_search_query, make_gather_query
from sourmash.index import LinearIndex


def test_make_jaccard_search_query():
    search_obj = make_jaccard_search_query(threshold=0)

    assert search_obj.score_fn == search_obj.score_jaccard
    assert not search_obj.require_scaled
    assert search_obj.threshold == 0


def test_make_jaccard_search_query_cont():
    search_obj = make_jaccard_search_query(do_containment=True,
                                           threshold=0)

    assert search_obj.score_fn == search_obj.score_containment
    assert search_obj.require_scaled
    assert search_obj.threshold == 0


def test_make_jaccard_search_query_max_cont():
    search_obj = make_jaccard_search_query(do_max_containment=True,
                                           threshold=0)

    assert search_obj.score_fn == search_obj.score_max_containment
    assert search_obj.require_scaled
    assert search_obj.threshold == 0


def test_make_jaccard_search_query_best_only():
    search_obj = make_jaccard_search_query(best_only=True)

    assert search_obj.score_fn == search_obj.score_jaccard
    assert not search_obj.require_scaled
    assert type(search_obj) == search.JaccardSearchBestOnly


def test_make_jaccard_search_query_no_threshold_none():
    search_obj = make_jaccard_search_query(threshold=None)

    assert search_obj.score_fn == search_obj.score_jaccard
    assert not search_obj.require_scaled
    assert search_obj.threshold == 0


def test_make_jaccard_search_query_cont_and_max_cont():
    with pytest.raises(TypeError) as exc:
        search_obj = make_jaccard_search_query(do_containment=True,
                                               do_max_containment=True)

    assert str(exc.value) == "'do_containment' and 'do_max_containment' cannot both be True"


def test_cont_requires_scaled():
    search_obj = make_jaccard_search_query(do_containment=True)
    assert search_obj.require_scaled
    
    mh = MinHash(n=500, ksize=31)
    with pytest.raises(TypeError) as exc:
        search_obj.check_is_compatible(SourmashSignature(mh))
    assert str(exc.value) == "this search requires a scaled signature"




def test_search_requires_flat():
    search_obj = make_jaccard_search_query()
    
    mh = MinHash(n=500, ksize=31, track_abundance=True)
    with pytest.raises(TypeError) as exc:
        search_obj.check_is_compatible(SourmashSignature(mh))
    assert str(exc.value) == "this search cannot be done with an abund signature"


def test_score_jaccard_similarity():
    search_obj = make_jaccard_search_query()

    assert search_obj.score_fn(None, 100, None, 200) == 0.5


def test_score_jaccard_containment():
    search_obj = make_jaccard_search_query(do_containment=True)

    assert search_obj.score_fn(100, 50, None, 0) == 0.5


def test_score_jaccard_containment_zero_query_size():
    search_obj = make_jaccard_search_query(do_containment=True)

    assert search_obj.score_fn(0, 100, None, None) == 0


def test_score_jaccard_max_containment_1():
    search_obj = make_jaccard_search_query(do_max_containment=True)

    assert search_obj.score_fn(150, 75, 100, None) == 0.75


def test_score_jaccard_max_containment_2():
    search_obj = make_jaccard_search_query(do_max_containment=True)

    assert search_obj.score_fn(100, 75, 150, None) == 0.75


def test_score_jaccard_max_containment_zero_query_size():
    search_obj = make_jaccard_search_query(do_containment=True)

    assert search_obj.score_fn(0, 100, None, None) == 0


def test_collect():
    search_obj = make_jaccard_search_query(threshold=0)
    search_obj.collect(1.0, None)
    assert search_obj.threshold == 0


def test_collect_best_only():
    search_obj = make_jaccard_search_query(threshold=0, best_only=True)
    search_obj.collect(1.0, None)
    assert search_obj.threshold == 1.0


def test_make_gather_query():
    # test basic make_gather_query call
    mh = MinHash(n=0, ksize=31, scaled=1000)

    for i in range(100):
        mh.add_hash(i)

    search_obj = make_gather_query(mh, 5e4)

    assert search_obj.score_fn == search_obj.score_containment
    assert search_obj.require_scaled
    assert search_obj.threshold == 0.5


def test_make_gather_query_no_threshold():
    # test basic make_gather_query call
    mh = MinHash(n=0, ksize=31, scaled=1000)

    for i in range(100):
        mh.add_hash(i)

    search_obj = make_gather_query(mh, None)

    assert search_obj.score_fn == search_obj.score_containment
    assert search_obj.require_scaled
    assert search_obj.threshold == 0


def test_make_gather_query_num_minhash():
    # will fail on non-scaled minhash
    mh = MinHash(n=500, ksize=31)

    for i in range(100):
        mh.add_hash(i)

    with pytest.raises(TypeError) as exc:
        search_obj = make_gather_query(mh, 5e4)

    assert str(exc.value) == "query signature must be calculated with scaled"


def test_make_gather_query_empty_minhash():
    # will fail on non-scaled minhash
    mh = MinHash(n=0, ksize=31, scaled=1000)

    for i in range(100):
        mh.add_hash(i)

    with pytest.raises(TypeError) as exc:
        search_obj = make_gather_query(mh, -1)

    assert str(exc.value) == "threshold_bp must be non-negative"


def test_make_gather_query_high_threshold():
    # will fail on non-scaled minhash
    mh = MinHash(n=0, ksize=31, scaled=1000)

    for i in range(100):
        mh.add_hash(i)

    # effective threshold > 1; raise ValueError
    with pytest.raises(ValueError):
        search_obj = make_gather_query(mh, 200000)


class FakeIndex(LinearIndex):
    _signatures = []
    filename = "something_or_other"

    def __init__(self, validator_fn):
        self.validator = validator_fn

    def find(self, search_fn, query, *args, **kwargs):
        if self.validator:
            self.validator(search_fn, query, args, kwargs)
        else:
            assert 0, "what are we even doing here?"
        return []


def test_index_search_passthru():
    # check that kwargs are passed through from 'search' to 'find'
    query = None

    def validate_kwarg_passthru(search_fn, query, args, kwargs):
        assert "this_kw_arg" in kwargs
        assert kwargs["this_kw_arg"] == 5

    idx = FakeIndex(validate_kwarg_passthru)

    idx.search(query, threshold=0.0, this_kw_arg=5)


def test_index_gather_passthru():
    # check that kwargs are passed through from 'gather' to 'find'
    query = None

    def validate_kwarg_passthru(search_fn, query, args, kwargs):
        assert "this_kw_arg" in kwargs
        assert kwargs["this_kw_arg"] == 5

    idx = FakeIndex(validate_kwarg_passthru)

    idx.search(query, threshold=0.0, this_kw_arg=5)
