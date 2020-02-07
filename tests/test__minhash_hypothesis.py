import pytest

from hypothesis import given, example
import hypothesis.strategies as st

from sourmash import MinHash
from sourmash._minhash import get_max_hash_for_scaled


@given(st.lists(st.integers(min_value=0, max_value=2**64 - 1), min_size=10, max_size=1000),
       st.lists(st.integers(min_value=0, max_value=2**64 - 1), min_size=10, max_size=1000),
       st.integers(min_value=10, max_value=1000))
@example([1, 2], [3, 4], 2)
def test_set_abundance_num_hypothesis(hashes, abundances, sketch_size):
    a = MinHash(sketch_size, 10, track_abundance=True)
    oracle = dict(zip(hashes, abundances))

    a.set_abundances(oracle)

    mins = a.get_mins(with_abundance=True)
    size = min(len(oracle), sketch_size)
    assert len(mins) == size

    for k, v in mins.items():
        assert oracle[k] == v


@given(st.lists(st.integers(min_value=0, max_value=2**64 - 1), min_size=10, max_size=1000),
       st.lists(st.integers(min_value=0, max_value=2**64 - 1), min_size=10, max_size=1000),
       st.integers(min_value=1000, max_value=10000))
def test_set_abundance_scaled_hypothesis(hashes, abundances, scaled):
    a = MinHash(0, 10, track_abundance=True, scaled=scaled)
    oracle = dict(zip(hashes, abundances))

    a.set_abundances(oracle)

    max_hash = get_max_hash_for_scaled(scaled)
    below_max_hash = sum(1 for k in oracle if k <= max_hash)

    mins = a.get_mins(with_abundance=True)
    assert len(mins) == below_max_hash

    for k, v in mins.items():
        assert oracle[k] == v
        assert k <= max_hash
