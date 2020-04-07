# This file is part of sourmash, https://github.com/dib-lab/sourmash/, and is
# Copyright (C) 2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: titus@idyll.org
# pylint: disable=missing-docstring,protected-access

from __future__ import print_function
from __future__ import absolute_import, unicode_literals

import itertools
import pickle
import math

import pytest

import sourmash
from sourmash._minhash import (
    MinHash,
    hash_murmur,
    get_scaled_for_max_hash,
    get_max_hash_for_scaled,
)
from sourmash import signature

from . import sourmash_tst_utils as utils

# add:
# * get default params from Python
# * keyword args for minhash constructor
# * trap error from handing protein/non-DNA to a DNA MH
# * fail on untagged/unloaded countgraph
# * nan on empty minhash
# * define equals


def test_basic_dna(track_abundance):
    # verify that MHs of size 1 stay size 1, & act properly as bottom sketches.
    mh = MinHash(1, 4, track_abundance=track_abundance)
    assert mh.moltype == 'DNA'

    mh.add_sequence('ATGC')
    a = mh.get_mins()

    mh.add_sequence('GCAT')             # this will not get added; hash > ATGC
    b = mh.get_mins()

    print(a, b)
    assert a == b
    assert len(b) == 1
    assert a[0] == b[0] == 12415348535738636339


def test_div_zero(track_abundance):
    # verify that empty MHs do not yield divide by zero errors for similarity
    mh = MinHash(1, 4, track_abundance=track_abundance)
    mh2 = mh.copy_and_clear()

    mh.add_sequence('ATGC')
    assert mh.similarity(mh2) == 0
    assert mh2.similarity(mh) == 0


def test_div_zero_contained(track_abundance):
    # verify that empty MHs do not yield divide by zero errors for contained_by
    mh = MinHash(1, 4, track_abundance=track_abundance)
    mh2 = mh.copy_and_clear()

    mh.add_sequence('ATGC')
    assert mh.contained_by(mh2) == 0
    assert mh2.contained_by(mh) == 0


def test_bytes_dna(track_abundance):
    mh = MinHash(1, 4, track_abundance=track_abundance)
    mh.add_sequence('ATGC')
    mh.add_sequence(b'ATGC')
    mh.add_sequence('ATGC')
    a = mh.get_mins()

    mh.add_sequence('GCAT')             # this will not get added; hash > ATGC
    mh.add_sequence(b'GCAT')             # this will not get added; hash > ATGC
    mh.add_sequence('GCAT')             # this will not get added; hash > ATGC
    b = mh.get_mins()

    print(a, b)
    assert a == b
    assert len(b) == 1


def test_bytes_protein_dayhoff(track_abundance, dayhoff):
    # verify that we can hash protein/aa sequences
    mh = MinHash(10, 6, True, dayhoff=dayhoff, hp=False,
                 track_abundance=track_abundance)

    expected_moltype = 'protein'
    if dayhoff:
        expected_moltype = 'dayhoff'
    assert mh.moltype == expected_moltype

    mh.add_protein('AGYYG')
    mh.add_protein('AGYYG')
    mh.add_protein(b'AGYYG')

    assert len(mh.get_mins()) == 4


def test_protein_dayhoff(track_abundance, dayhoff):
    # verify that we can hash protein/aa sequences
    mh = MinHash(10, 6, True, dayhoff=dayhoff, hp=False, track_abundance=track_abundance)
    mh.add_protein('AGYYG')

    assert len(mh.get_mins()) == 4


def test_bytes_protein_hp(track_abundance, hp):
    # verify that we can hash protein/aa sequences
    mh = MinHash(10, 6, True, dayhoff=False, hp=hp, track_abundance=track_abundance)
    expected_moltype = 'protein'
    if hp:
        expected_moltype = 'hp'
    assert mh.moltype == expected_moltype

    mh.add_protein('AGYYG')
    mh.add_protein(u'AGYYG')
    mh.add_protein(b'AGYYG')

    if hp:
        assert len(mh.get_mins()) == 1
    else:
        assert len(mh.get_mins()) == 4


def test_protein_hp(track_abundance, hp):
    # verify that we can hash protein/aa sequences
    mh = MinHash(10, 6, True, dayhoff=False, hp=hp, track_abundance=track_abundance)
    mh.add_protein('AGYYG')

    if hp:
        assert len(mh.get_mins()) == 1
    else:
        assert len(mh.get_mins()) == 4


def test_translate_codon(track_abundance):
    # Ensure that translation occurs properly
    mh = MinHash(10, 6, is_protein=True)
    assert mh.moltype == 'protein'

    assert "S" == mh.translate_codon('TCT')
    assert "S" == mh.translate_codon('TC')
    assert "X" == mh.translate_codon("T")

    with pytest.raises(ValueError):
        mh.translate_codon("")
        mh.translate_codon("TCTA")


def test_dayhoff(track_abundance):
    # verify that we can hash to dayhoff-encoded protein/aa sequences
    mh_dayhoff = MinHash(10, 6, is_protein=True,
                         dayhoff=True, hp=False, track_abundance=track_abundance)
    mh_dayhoff.add_sequence('ACTGAC')

    assert len(mh_dayhoff.get_mins()) == 2
    # verify that dayhoff-encoded hashes are different from protein/aa hashes
    mh_protein = MinHash(10, 6, is_protein=True, track_abundance=track_abundance)
    mh_protein.add_sequence('ACTGAC')

    assert len(mh_protein.get_mins()) == 2
    assert mh_protein.get_mins() != mh_dayhoff.get_mins()


def test_hp(track_abundance):
    # verify that we can hash to hp-encoded protein/aa sequences
    mh_hp = MinHash(10, 6, is_protein=True,
                    dayhoff=False, hp=True, track_abundance=track_abundance)
    assert mh_hp.moltype == 'hp'

    mh_hp.add_sequence('ACTGAC')

    assert len(mh_hp.get_mins()) == 2
    # verify that hp-encoded hashes are different from protein/aa hashes
    mh_protein = MinHash(10, 6, is_protein=True, track_abundance=track_abundance)
    mh_protein.add_sequence('ACTGAC')

    assert len(mh_protein.get_mins()) == 2
    assert mh_protein.get_mins() != mh_hp.get_mins()


def test_protein_short(track_abundance):
    # verify that we can hash protein/aa sequences
    mh = MinHash(10, 9, True, track_abundance=track_abundance)
    mh.add_protein('AG')

    assert len(mh.get_mins()) == 0, mh.get_mins()


def test_size_limit(track_abundance):
    # test behavior with size limit of 3
    mh = MinHash(3, 4, track_abundance=track_abundance)
    mh.add_hash(10)
    mh.add_hash(20)
    mh.add_hash(30)
    assert mh.get_mins() == [10, 20, 30]
    mh.add_hash(5) # -> should push 30 off end
    assert mh.get_mins() == [5, 10, 20]


def test_max_hash(track_abundance):
    # test behavior with max_hash
    mh = MinHash(0, 4, track_abundance=track_abundance, max_hash=35)
    mh.add_hash(10)
    mh.add_hash(20)
    mh.add_hash(30)
    assert mh.get_mins() == [10, 20, 30]
    mh.add_hash(40)
    assert mh.get_mins() == [10, 20, 30]
    mh.add_hash(36)
    assert mh.get_mins() == [10, 20, 30]


def test_scaled(track_abundance):
    # test behavior with scaled (alt to max_hash)
    scaled = get_scaled_for_max_hash(35)
    print('XX', scaled, get_max_hash_for_scaled(scaled))
    mh = MinHash(0, 4, track_abundance=track_abundance, scaled=scaled)
    assert mh.max_hash == 35

    mh.add_hash(10)
    mh.add_hash(20)
    mh.add_hash(30)
    assert mh.get_mins() == [10, 20, 30]
    mh.add_hash(40)
    assert mh.get_mins() == [10, 20, 30]
    mh.add_hash(36)
    assert mh.get_mins() == [10, 20, 30]


def test_no_scaled(track_abundance):
    # no 'scaled', num=0 - should fail
    with pytest.raises(ValueError):
        mh = MinHash(0, 4, track_abundance=track_abundance)


def test_max_hash_conversion():
    SCALED=100000
    max_hash = get_max_hash_for_scaled(SCALED)
    new_scaled = get_scaled_for_max_hash(max_hash)
    assert new_scaled == SCALED


def test_max_hash_and_scaled_zero():
    max_hash = get_max_hash_for_scaled(0)
    new_scaled = get_scaled_for_max_hash(0)
    assert max_hash == new_scaled
    assert max_hash == 0


def test_max_hash_and_scaled_error(track_abundance):
    # test behavior when supplying both max_hash and scaled
    with pytest.raises(ValueError):
        mh = MinHash(0, 4, track_abundance=track_abundance, max_hash=35,
                     scaled=5)


def test_max_hash_cannot_limit(track_abundance):
    # make sure you can't set both max_n and max_hash.
    with pytest.raises(ValueError):
        mh = MinHash(2, 4, track_abundance=track_abundance, max_hash=35)


def test_no_downsample_scaled_if_n(track_abundance):
    # make sure you can't set max_n and then downsample scaled
    mh = MinHash(2, 4, track_abundance=track_abundance)
    with pytest.raises(ValueError) as excinfo:
        mh.downsample_scaled(100000000)

    assert 'cannot downsample a standard MinHash' in str(excinfo.value)


def test_scaled(track_abundance):
    # make sure you can't set both max_n and scaled.
    with pytest.raises(ValueError):
        mh = MinHash(2, 4, track_abundance=track_abundance, scaled=2)


def test_mh_jaccard_similarity():
    # check actual Jaccard value for a non-trivial case
    a = MinHash(0, 20, max_hash=50, track_abundance=False)
    b = MinHash(0, 20, max_hash=50, track_abundance=False)
    a.add_many([1, 3, 5, 8])
    b.add_many([1, 3, 5, 6, 8, 10])

    assert a.similarity(b) == 4. / 6.


def test_mh_similarity_downsample_jaccard_value():
    # check jaccard value after downsampling

    # max_hash = 50
    a = MinHash(0, 20, max_hash=50, track_abundance=False)
    # max_hash = 100
    b = MinHash(0, 20, max_hash=100, track_abundance=False)

    a.add_many([1, 3, 5, 8, 70])
    b.add_many([1, 3, 5, 6, 8, 10, 70 ])

    # the hash=70 will be truncated by downsampling
    assert a.similarity(b, downsample=True) == 4. / 6.


def test_mh_angular_similarity():
    # check actual angular similarity for a non-trivial case, taken from:
    # https://www.sciencedirect.com/topics/computer-science/cosine-similarity
    # note: angular similarity is 1 - 2*(acos(sim) / pi), when elements
    # are always positive (https://en.wikipedia.org/wiki/Cosine_similarity)
    a = MinHash(0, 20, max_hash=50, track_abundance=True)
    b = MinHash(0, 20, max_hash=50, track_abundance=True)
    a.set_abundances({ 1:5, 3:3, 5:2, 8:2})
    b.set_abundances({ 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 })

    cos_sim = 0.9356
    angular_sim = 1 - 2*math.acos(cos_sim) / math.pi
    assert round(angular_sim, 4) == 0.7703

    assert round(a.similarity(b), 4) == round(angular_sim, 4)


def test_mh_angular_similarity_2():
    # check actual angular similarity for a second non-trivial case
    a = MinHash(0, 20, max_hash=100, track_abundance=True)
    b = MinHash(0, 20, max_hash=100, track_abundance=True)
    a.set_abundances({ 1:5, 3:3, 5:2, 8:2, 70:70 })
    b.set_abundances({ 1:3, 3:2, 5:1, 6:1, 8:1, 10:1, 70:70 })

    assert round(a.similarity(b), 4) == 0.9728

    # ignore_abundance => jaccard
    assert a.similarity(b, ignore_abundance=True) == 5. / 7.


def test_mh_similarity_downsample_angular_value():
    # test downsample=True argument to MinHash.similarity

    # max_hash = 50
    a = MinHash(0, 20, max_hash=50, track_abundance=True)
    # max_hash = 100
    b = MinHash(0, 20, max_hash=100, track_abundance=True)

    a.set_abundances({ 1:5, 3:3, 5:2, 8:2, 70:70 })
    b.set_abundances({ 1:3, 3:2, 5:1, 6:1, 8:1, 10:1, 70:70 })

    # the hash=70 will be truncated by downsampling
    sim = a.similarity(b, downsample=True)
    assert round(sim, 4) == 0.7703

    # with ignore_abundance, will be equal to jaccard
    jaccard = a.similarity(b, downsample=True, ignore_abundance=True)
    assert jaccard == 4. / 6.


def test_mh_similarity_downsample_true(track_abundance):
    # verify sim(a, b) == sim(b, a), with and without ignore_abundance

    # max_hash = 50
    a = MinHash(0, 20, max_hash=50, track_abundance=track_abundance)
    # max_hash = 100
    b = MinHash(0, 20, max_hash=100, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }
    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    # downsample=True => no error; values should match either way
    x = a.similarity(b, ignore_abundance=True, downsample=True)
    y = b.similarity(a, ignore_abundance=True, downsample=True)
    assert x == y

    # downsample=True => no error; values should match either way
    x = a.similarity(b, ignore_abundance=False, downsample=True)
    y = b.similarity(a, ignore_abundance=False, downsample=True)
    assert x == y


def test_mh_similarity_downsample_errors(track_abundance):
    # test downsample=False (default) argument to MinHash.similarity

    # max_hash = 50
    a = MinHash(0, 20, max_hash=50, track_abundance=track_abundance)
    # max_hash = 100
    b = MinHash(0, 20, max_hash=100, track_abundance=track_abundance)

    a_values = { 1:5, 3:3, 5:2, 8:2}
    b_values = { 1:3, 3:2, 5:1, 6:1, 8:1, 10:1 }
    if track_abundance:
        a.set_abundances(a_values)
        b.set_abundances(b_values)
    else:
        a.add_many(a_values.keys())
        b.add_many(b_values.keys())

    # error, incompatible max hash
    with pytest.raises(ValueError) as e:
        a.similarity(b, ignore_abundance=True)   # downsample=False
    assert 'mismatch in scaled; comparison fail' in str(e.value)

    with pytest.raises(ValueError) as e:
        a.similarity(b, ignore_abundance=False)  # downsample=False
    assert 'mismatch in scaled; comparison fail' in str(e.value)

    with pytest.raises(ValueError) as e:
        b.similarity(a, ignore_abundance=True)   # downsample=False
    assert 'mismatch in scaled; comparison fail' in str(e.value)

    with pytest.raises(ValueError) as e:
        b.similarity(a, ignore_abundance=False)  # downsample=false
    assert 'mismatch in scaled; comparison fail' in str(e.value)


def test_basic_dna_bad(track_abundance):
    # test behavior on bad DNA
    mh = MinHash(1, 4, track_abundance=track_abundance)

    with pytest.raises(ValueError) as e:
        mh.add_sequence('ATGR')
    print(e)

    assert 'invalid DNA character in input k-mer: ATGR' in str(e.value)


def test_basic_dna_bad_2(track_abundance):
    # test behavior on bad DNA
    mh = MinHash(1, 6, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        mh.add_protein('YYYY')


def test_basic_dna_bad_force(track_abundance):
    # test behavior on bad DNA; use 100 so multiple hashes get added.
    mh = MinHash(100, 4, track_abundance=track_abundance)
    assert len(mh.get_mins()) == 0
    mh.add_sequence('ATGN', True)     # ambiguous kmer skipped.
    assert len(mh.get_mins()) == 0
    mh.add_sequence('AATGN', True)    # but good k-mers still used.
    assert len(mh.get_mins()) == 1
    mh.add_sequence('AATG', True)     # checking that right kmer was added
    assert len(mh.get_mins()) == 1    # (only 1 hash <- this is a dup)


def test_basic_dna_bad_force_2(track_abundance):
    # test behavior on bad DNA
    mh = MinHash(100, 4, track_abundance=track_abundance)
    assert len(mh.get_mins()) == 0
    mh.add_sequence('AAGNCGG', True)     # ambiguous kmers skipped.
    assert len(mh.get_mins()) == 0
    mh.add_sequence('AATGNGCGG', True)  # ambiguous kmers skipped.
    assert len(mh.get_mins()) == 2
    mh.add_sequence('AATG', True)        # checking that right kmers were added
    mh.add_sequence('GCGG', True)
    assert len(mh.get_mins()) == 2       # (only 2 hashes should be there)


def test_consume_lowercase(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    b = MinHash(20, 10, track_abundance=track_abundance)

    a.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA'.lower())
    b.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')

    assert a.compare(b) == 1.0
    assert b.compare(b) == 1.0
    assert b.compare(a) == 1.0
    assert a.compare(a) == 1.0


def test_compare_1(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    b = MinHash(20, 10, track_abundance=track_abundance)

    a.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    b.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')

    assert a.compare(b) == 1.0
    assert b.compare(b) == 1.0
    assert b.compare(a) == 1.0
    assert a.compare(a) == 1.0

    # add same sequence again
    b.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    assert a.compare(b) == 1.0
    assert b.compare(b) == 1.0
    assert b.compare(a) == 1.0
    assert a.compare(a) == 1.0


    b.add_sequence('GATTGGTGCACACTTAACTGGGTGCCGCGCTGGTGCTGATCCATGAAGTT')
    x = a.compare(b)
    assert x >= 0.3, x

    x = b.compare(a)
    assert x >= 0.3, x
    assert a.compare(a) == 1.0
    assert b.compare(b) == 1.0


def test_intersection_errors(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    b = MinHash(20, 10, track_abundance=track_abundance)
    c = MinHash(30, 10, track_abundance=track_abundance)

    a.add_sequence("TGCCGCCCAGCA")
    b.add_sequence("TGCCGCCCAGCA")

    common = set(a.get_mins())
    combined_size = 3

    intersection, size = a.intersection(b, in_common=False)
    assert intersection == set()
    assert combined_size == size

    with pytest.raises(TypeError):
        a.intersection(set())

    with pytest.raises(TypeError):
        a.intersection(c)


# this filter doesn't work, but leaving it in pour encourages les autres.
@pytest.mark.filterwarnings("ignore")
def test_intersection_1(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    b = MinHash(20, 10, track_abundance=track_abundance)

    a.add_sequence('TGCCGCCCAGCA')
    b.add_sequence('TGCCGCCCAGCA')

    common = set(a.get_mins())
    combined_size = 3

    intersection, size = a.intersection(b, in_common=True)
    assert intersection == common
    assert combined_size == size

    intersection, size = b.intersection(b, in_common=True)
    assert intersection == common
    assert combined_size == size

    intersection, size = b.intersection(a, in_common=True)
    assert intersection == common
    assert combined_size == size

    intersection, size = a.intersection(a, in_common=True)
    assert intersection == common
    assert combined_size == size

    # add same sequence again
    b.add_sequence('TGCCGCCCAGCA')

    intersection, size = a.intersection(b, in_common=True)
    assert intersection == common
    assert combined_size == size

    intersection, size = b.intersection(b, in_common=True)
    assert intersection == common
    assert combined_size == size

    intersection, size = b.intersection(a, in_common=True)
    assert intersection == common
    assert combined_size == size

    intersection, size = a.intersection(a, in_common=True)
    assert intersection == common
    assert combined_size == size

    a.add_sequence('GTCCGCCCAGTGA')
    b.add_sequence('GTCCGCCCAGTGG')

    new_in_common = set(a.get_mins()).intersection(set(b.get_mins()))
    new_combined_size = 8

    intersection, size = a.intersection(b, in_common=True)
    assert intersection == new_in_common
    assert size == new_combined_size

    intersection, size = b.intersection(a, in_common=True)
    assert intersection == new_in_common
    assert size == new_combined_size

    intersection, size = a.intersection(a, in_common=True)
    assert intersection == set(a.get_mins())

    intersection, size = b.intersection(b, in_common=True)
    assert intersection == set(b.get_mins())


def test_mh_copy(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)

    a.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    b = a.__copy__()
    assert b.compare(a) == 1.0


def test_mh_len(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)

    assert len(a) == 20
    a.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    assert len(a) == 20


def test_mh_len(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    assert a.get_mins() == list(range(0, 40, 2))


def test_mh_unsigned_long_long(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    a.add_hash(9227159859419181011)        # too big for a C long int.
    assert 9227159859419181011 in a.get_mins()


def test_mh_count_common(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    assert a.count_common(b) == 10
    assert b.count_common(a) == 10


def test_mh_count_common_diff_protein(track_abundance):
    a = MinHash(20, 5, False, track_abundance=track_abundance)
    b = MinHash(20, 5, True, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a.count_common(b)


def test_mh_count_common_diff_maxhash(track_abundance):
    a = MinHash(0, 5, False, track_abundance=track_abundance, max_hash=1)
    b = MinHash(0, 5, True, track_abundance=track_abundance, max_hash=2)

    with pytest.raises(ValueError):
        a.count_common(b)


def test_mh_count_common_diff_seed(track_abundance):
    a = MinHash(20, 5, False, track_abundance=track_abundance, seed=1)
    b = MinHash(20, 5, True, track_abundance=track_abundance, seed=2)

    with pytest.raises(ValueError):
        a.count_common(b)


def test_mh_count_common_diff_ksize(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance)
    b = MinHash(20, 6, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a.count_common(b)


def test_mh_count_common_notmh(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance)
    b = set()

    with pytest.raises(TypeError):
        a.count_common(b)


def test_mh_downsample_n_error(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    with pytest.raises(ValueError):
        a.downsample_n(30)


def test_mh_jaccard_asymmetric_num(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    # different size: 10
    b = MinHash(10, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    assert a.count_common(b) == 10
    assert b.count_common(a) == 10

    with pytest.raises(TypeError):
        a.compare(b)

    a = a.downsample_n(10)
    assert a.compare(b) == 0.5
    assert b.compare(a) == 0.5


def test_mh_merge_typeerror(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    with pytest.raises(TypeError):
        a.merge(set())


def test_mh_merge(track_abundance):
    # test merging two identically configured minhashes
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.merge(b)
    d = b.merge(a)

    assert len(c) == len(d)
    assert c.get_mins() == d.get_mins()
    assert c.compare(d) == 1.0
    assert d.compare(c) == 1.0


def test_mh_merge_empty_num(track_abundance):
    # test merging two identically configured minhashes, one empty
    a = MinHash(20, 10, track_abundance=track_abundance)

    b = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.merge(b)
    d = b.merge(a)

    assert len(c)
    assert len(c) == len(d)
    assert c.get_mins() == d.get_mins()
    assert c.compare(d) == 1.0
    assert d.compare(c) == 1.0


def test_mh_merge_empty_scaled(track_abundance):
    # test merging two identically configured minhashes, one empty
    a = MinHash(0, 10, scaled=1, track_abundance=track_abundance)

    b = MinHash(0, 10, scaled=1, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.merge(b)
    d = b.merge(a)

    assert len(c)
    assert len(c) == len(d)
    assert c.get_mins() == d.get_mins()
    assert c.compare(d) == 1.0
    assert d.compare(c) == 1.0


def test_mh_merge_check_length(track_abundance):
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.merge(b)
    assert len(c.get_mins()) == 20


def test_mh_merge_check_length2(track_abundance):
    # merged MH doesn't have full number of elements
    a = MinHash(4, 10, track_abundance=track_abundance)
    a.add_hash(3)
    a.add_hash(1)
    a.add_hash(4)

    b = MinHash(4, 10, track_abundance=track_abundance)
    b.add_hash(3)
    b.add_hash(1)
    b.add_hash(4)

    c = a.merge(b)
    assert len(c.get_mins()) == 3

def test_mh_asymmetric_merge(track_abundance):
    # test merging two asymmetric (different size) MHs
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    # different size: 10
    b = MinHash(10, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.merge(b)
    d = b.merge(a)

    assert len(a) == 20
    assert len(b) == 10
    assert len(c) == len(a)
    assert len(d) == len(b)

    # can't compare different sizes without downsampling
    with pytest.raises(TypeError):
        d.compare(a)

    a = a.downsample_n(d.num)
    print(a.get_mins())
    print(d.get_mins())
    assert d.compare(a) == 1.0

    c = c.downsample_n(b.num)
    assert c.compare(b) == 1.0


def test_mh_inplace_concat_asymmetric(track_abundance):
    # test merging two asymmetric (different size) MHs
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    # different size: 10
    b = MinHash(10, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.__copy__()
    c += b

    d = b.__copy__()
    d += a

    assert len(a) == 20
    assert len(b) == 10
    assert len(c) == len(a)
    assert len(d) == len(b)

    try:
        d.compare(a)
    except TypeError as exc:
        assert 'must have same num' in str(exc)

    a = a.downsample_n(d.num)
    assert d.compare(a) == 1.0 # see: d += a, above.

    c = c.downsample_n(b.num)
    assert c.compare(b) == 0.5


def test_mh_inplace_concat(track_abundance):
    # test merging two identically configured minhashes
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    c = a.__copy__()
    c += b
    d = b.__copy__()
    d += a

    assert len(c) == len(d)
    assert c.get_mins() == d.get_mins()
    assert c.compare(d) == 1.0
    assert d.compare(c) == 1.0


def test_mh_merge_diff_protein(track_abundance):
    a = MinHash(20, 5, False, track_abundance=track_abundance)
    b = MinHash(20, 5, True, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a.merge(b)


def test_mh_merge_diff_ksize(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance)
    b = MinHash(20, 6, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a.merge(b)


def test_mh_compare_diff_protein(track_abundance):
    a = MinHash(20, 5, False, track_abundance=track_abundance)
    b = MinHash(20, 5, True, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a.compare(b)


def test_mh_compare_diff_ksize(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance)
    b = MinHash(20, 6, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a.compare(b)


def test_mh_compare_diff_seed(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance, seed=1)
    b = MinHash(20, 5, track_abundance=track_abundance, seed=2)

    with pytest.raises(ValueError):
        a.compare(b)


def test_mh_compare_diff_max_hash(track_abundance):
    a = MinHash(0, 5, track_abundance=track_abundance, max_hash=5)
    b = MinHash(0, 5, track_abundance=track_abundance, max_hash=10)

    with pytest.raises(ValueError):
        a.compare(b)


def test_mh_concat_diff_protein(track_abundance):
    a = MinHash(20, 5, False, track_abundance=track_abundance)
    b = MinHash(20, 5, True, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a += b


def test_mh_concat_diff_ksize(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance)
    b = MinHash(20, 6, track_abundance=track_abundance)

    with pytest.raises(ValueError):
        a += b


def test_mh_concat_diff_max_hash(track_abundance):
    a = MinHash(0, 5, track_abundance=track_abundance, max_hash=5)
    b = MinHash(0, 5, track_abundance=track_abundance, max_hash=10)

    with pytest.raises(ValueError):
        a += b


def test_mh_concat_diff_seed(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance, seed=1)
    b = MinHash(20, 5, track_abundance=track_abundance, seed=2)

    with pytest.raises(ValueError):
        a += b


def test_short_sequence(track_abundance):
    a = MinHash(20, 5, track_abundance=track_abundance)
    a.add_sequence('GGGG')
    # adding a short sequence should fail silently
    assert len(a.get_mins()) == 0


def test_bytes_murmur():
    x = hash_murmur("ACG")
    assert x == 1731421407650554201

    x = hash_murmur(b"ACG")
    assert x == 1731421407650554201

    x = hash_murmur(u"ACG")
    assert x == 1731421407650554201


def test_murmur():
    x = hash_murmur("ACG")
    assert x == 1731421407650554201

    try:
        x = hash_murmur()
        assert 0, "hash_murmur requires an argument"
    except TypeError:
        pass

    x = hash_murmur("ACG", 42)
    assert x == 1731421407650554201

    y = hash_murmur("ACG", 43)
    assert y != x


def test_abundance_simple():
    a = MinHash(20, 5, False, track_abundance=True)

    a.add_sequence('AAAAA')
    assert a.get_mins() == [2110480117637990133]
    assert a.get_mins(with_abundance=True) == {2110480117637990133: 1}

    a.add_sequence('AAAAA')
    assert a.get_mins() == [2110480117637990133]
    assert a.get_mins(with_abundance=True) == {2110480117637990133: 2}


def test_abundance_simple_2():
    a = MinHash(20, 5, False, track_abundance=True)
    b = MinHash(20, 5, False, track_abundance=True)

    a.add_sequence('AAAAA')
    assert a.get_mins() == [2110480117637990133]
    assert a.get_mins(with_abundance=True) == {2110480117637990133: 1}

    a.add_sequence('AAAAA')
    assert a.get_mins() == [2110480117637990133]
    assert a.get_mins(with_abundance=True) == {2110480117637990133: 2}

    b.add_sequence('AAAAA')
    assert a.count_common(b) == 1


def test_abundance_count_common():
    a = MinHash(20, 5, False, track_abundance=True)
    b = MinHash(20, 5, False, track_abundance=False)

    a.add_sequence('AAAAA')
    a.add_sequence('AAAAA')
    assert a.get_mins() == [2110480117637990133]
    assert a.get_mins(with_abundance=True) == {2110480117637990133: 2}

    b.add_sequence('AAAAA')
    b.add_sequence('GGGGG')
    assert a.count_common(b) == 1
    assert a.count_common(b) == b.count_common(a)

    assert b.get_mins(with_abundance=True) == [2110480117637990133,
                                               10798773792509008305]


def test_abundance_compare():
    a = MinHash(20, 10, track_abundance=True)
    b = MinHash(20, 10, track_abundance=False)

    a.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    b.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')

    assert a.compare(b) == 1.0
    assert b.compare(b) == 1.0
    assert b.compare(a) == 1.0
    assert a.compare(a) == 1.0

    # add same sequence again
    b.add_sequence('TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA')
    assert a.compare(b) == 1.0
    assert b.compare(b) == 1.0
    assert b.compare(a) == 1.0
    assert a.compare(a) == 1.0

    b.add_sequence('GATTGGTGCACACTTAACTGGGTGCCGCGCTGGTGCTGATCCATGAAGTT')
    x = a.compare(b)
    assert x >= 0.3, x

    x = b.compare(a)
    assert x >= 0.3, x
    assert a.compare(a) == 1.0
    assert b.compare(b) == 1.0


def test_set_abundance():
    a = MinHash(20, 10, track_abundance=False)

    with pytest.raises(RuntimeError) as e:
        a.set_abundances({1: 3, 2: 4})

    assert "track_abundance=True when constructing" in e.value.args[0]


def test_set_abundance_2():
    datapath = utils.get_test_data("genome-s12.fa.gz.sig")
    sig = sourmash.load_one_signature(datapath,
                                      ksize=30,
                                      select_moltype='dna')
    new_mh = sig.minhash.copy_and_clear()
    mins = sig.minhash.get_mins()
    mins = {k: 1 for k in mins}
    new_mh.track_abundance = True
    new_mh.set_abundances(mins)

    assert new_mh.get_mins(with_abundance=True) == mins


def test_reset_abundance_initialized():
    a = MinHash(1, 4, track_abundance=True)
    a.add_sequence('ATGC')

    # If we had a minhash with abundances and drop it, this shouldn't fail.
    # Convert from Abundance to Regular MinHash
    a.track_abundance = False

    assert a.get_mins(with_abundance=True) == [12415348535738636339]


def test_set_abundance_initialized():
    a = MinHash(1, 4, track_abundance=False)
    a.add_sequence('ATGC')

    with pytest.raises(RuntimeError) as e:
        a.track_abundance = True

    assert "Can only set track_abundance=True if the MinHash is empty" in e.value.args[0]


def test_reviving_minhash():
    # simulate reading a MinHash from disk
    mh = MinHash(0, 21, max_hash=184467440737095520, seed=42,
                 track_abundance=False)
    mins = (28945103950853965, 74690756200987412, 82962372765557409,
            93503551367950366, 106923350319729608, 135116761470196737,
            160165359281648267, 162390811417732001, 177939655451276972)

    for m in mins:
        mh.add_hash(m)


def test_set_abundance_num():
    a = MinHash(2, 10, track_abundance=True)

    a.set_abundances({1: 3, 2: 4})

    assert a.get_mins(with_abundance=True) == {1: 3, 2: 4}


def test_mh_copy_and_clear(track_abundance):
    # test basic creation of new, empty MinHash
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = a.copy_and_clear()
    assert a.ksize == b.ksize
    assert b.num == a.num
    assert b.max_hash == 0
    assert not b.is_protein
    assert b.track_abundance == track_abundance
    assert b.seed == a.seed
    assert len(b.get_mins()) == 0
    assert a.scaled == b.scaled
    assert b.scaled == 0


def test_mh_copy_and_clear_with_max_hash(track_abundance):
    # test basic creation of new, empty MinHash w/max_hash param set
    a = MinHash(0, 10, track_abundance=track_abundance, max_hash=20)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = a.copy_and_clear()
    assert a.ksize == b.ksize
    assert b.num == a.num
    assert b.max_hash == 20
    assert not b.is_protein
    assert b.track_abundance == track_abundance
    assert b.seed == a.seed
    assert len(b.get_mins()) == 0
    assert a.scaled == b.scaled
    assert b.scaled != 0


def test_scaled_property(track_abundance):
    scaled = 10000
    a = MinHash(0, 10, track_abundance=track_abundance,
                max_hash=round(2**64 / scaled))
    assert a.scaled == scaled


def test_mh_subtract(track_abundance):
    # test subtracting two identically configured minhashes
    a = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = MinHash(20, 10, track_abundance=track_abundance)
    for i in range(0, 80, 4):
        b.add_hash(i)

    assert a.subtract_mins(b) == set(range(2, 40, 4))


def test_pickle_max_hash(track_abundance):
    a = MinHash(0, 10, track_abundance=track_abundance, max_hash=20)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = pickle.loads(pickle.dumps(a))
    assert a.ksize == b.ksize
    assert b.num == a.num
    assert b.max_hash == a.max_hash
    assert b.max_hash == 20
    assert not b.is_protein
    assert b.track_abundance == track_abundance
    assert b.seed == a.seed
    assert len(b.get_mins()) == len(a.get_mins())
    assert len(b.get_mins()) == 11
    assert a.scaled == b.scaled
    assert b.scaled != 0


def test_pickle_scaled(track_abundance):
    a = MinHash(0, 10, track_abundance=track_abundance, scaled=922337203685477632)
    for i in range(0, 40, 2):
        a.add_hash(i)

    b = pickle.loads(pickle.dumps(a))
    assert a.ksize == b.ksize
    assert b.num == a.num
    assert b.max_hash == a.max_hash
    assert b.max_hash == 20
    assert not b.is_protein
    assert b.track_abundance == track_abundance
    assert b.seed == a.seed
    assert len(b.get_mins()) == len(a.get_mins())
    assert len(b.get_mins()) == 11
    assert a.scaled == b.scaled
    assert b.scaled != 0


def test_minhash_abund_add():
    # this targets part of bug #319, a segfault caused by invalidation of
    # std::vector iterators upon vector resizing - in this case, there
    # was also a bug in inserting into the middle of mins when scaled was set.

    a = MinHash(0, 10, track_abundance=True, max_hash=5000)

    n = 0
    for i in range(10, 0, -1):
        a.add_hash(i)
        n += 1
        assert len(a.get_mins()) == n
        print(len(a.get_mins()))


def test_minhash_abund_capacity_increase():
    # this targets bug #319, a segfault caused by invalidation of
    # std::vector iterators upon vector resizing.

    # this should set capacity to 1000 - see KmerMinHash constructor call
    # to 'reserve' when n > 0 for specific parameter.
    a = MinHash(0, 10, track_abundance=True, max_hash=5000)

    # 1001 is dependent on the value passed to reserve (currently 1000).
    for i in range(1001, 0, -1):
        a.add_hash(i)


def test_minhash_abund_merge_flat():
    # this targets a segfault caused by trying to compute similarity
    # of a signature with abundance and a signature without abundance.
    # the correct behavior for now is to calculate simple Jaccard,
    # i.e. 'flatten' both of them.
    a = MinHash(0, 10, track_abundance=True, max_hash=5000)
    b = MinHash(0, 10, max_hash=5000)

    for i in range(0, 10, 2):
        a.add_hash(i)

    for j in range(0, 10, 3):
        b.add_hash(i)

    # these crashed, previously.
    assert a.similarity(b) == 0.2
    assert b.similarity(a) == 0.2


def test_minhash_abund_merge_flat_2():
    # this targets a segfault caused by trying to merge
    # a signature with abundance and a signature without abundance.

    a = MinHash(0, 10, track_abundance=True, max_hash=5000)
    b = MinHash(0, 10, max_hash=5000)

    for i in range(0, 10, 2):
        a.add_hash(i)

    for j in range(0, 10, 3):
        b.add_hash(i)

    a.merge(b)


def test_distance_matrix(track_abundance):
    import numpy

    siglist = [next(signature.load_signatures(utils.get_test_data(f)))
               for f in utils.SIG_FILES]

    D1 = numpy.zeros([len(siglist), len(siglist)])
    D2 = numpy.zeros([len(siglist), len(siglist)])

    for i, E in enumerate(siglist):
        for j, E2 in enumerate(siglist):
            if i < j:
                continue
            similarity = E.similarity(E2, track_abundance)
            D2[i][j] = similarity
            D2[j][i] = similarity

    for i, E in enumerate(siglist):
        for j, E2 in enumerate(siglist):
            D1[i][j] = E.similarity(E2, track_abundance)

    assert numpy.array_equal(D1, D2)


def test_remove_many(track_abundance):
    a = MinHash(0, 10, track_abundance=track_abundance, max_hash=5000)

    a.add_many(list(range(0, 100, 2)))

    orig_sig = signature.SourmashSignature(a)
    orig_md5 = orig_sig.md5sum()

    a.remove_many(list(range(0, 100, 3)))
    new_sig = signature.SourmashSignature(a)
    new_md5 = new_sig.md5sum()

    assert orig_md5 == "f1cc295157374f5c07cfca5f867188a1"
    assert new_md5 == "dd93fa319ef57f4a019c59ee1a8c73e2"
    assert orig_md5 != new_md5

    assert len(a) == 33
    assert all(c % 6 != 0 for c in a.get_mins())


def test_add_many(track_abundance):
    a = MinHash(0, 10, track_abundance=track_abundance, max_hash=5000)
    b = MinHash(0, 10, track_abundance=track_abundance, max_hash=5000)

    a.add_many(list(range(0, 100, 2)))
    a.add_many(list(range(0, 100, 2)))

    assert len(a) == 50
    assert all(c % 2 == 0 for c in a.get_mins())

    for h in range(0, 100, 2):
        b.add_hash(h)
        b.add_hash(h)

    assert len(b) == 50
    assert a == b


def test_set_abundances_huge():
    max_hash = 4000000
    a = MinHash(0, 10, track_abundance=True, max_hash=max_hash)

    hashes = list(range(max_hash))
    abundances = itertools.repeat(2)

    a.set_abundances(dict(zip(hashes, abundances)))
