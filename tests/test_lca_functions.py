# CTB: move other functions from test_lca.py into here.
from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import *


class FakeLCA_Database(object):
    def __init__(self):
        self._assignments = {}

    def _set_lineage_assignment(self, hashval, assignment):
        self._assignments[hashval] = assignment

    def get_lineage_assignments(self, hashval):
        if hashval in self._assignments:
            return self._assignments[hashval]
        else:
            return None


def test_gather_assignments_1():
    # test basic mechanics of gather_assignments function
    hashval = 12345678
    lin = lca_utils.make_lineage('a;b;c')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin ]))

    assignments = lca_utils.gather_assignments([hashval], [db])
    print(assignments)

    assert assignments[hashval] == set([ lin ])


def test_gather_assignments_2():
    # test basic mechanics of gather_assignments function with two lineages
    hashval = 12345678
    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))

    assignments = lca_utils.gather_assignments([hashval], [db])
    print(assignments)

    assert assignments[hashval] == set([ lin, lin2 ])


def test_gather_assignments_3():
    # test basic mechanics of gather_assignments function with two lineages
    # and two hashvals
    hashval = 12345678
    hashval2 = 87654321
    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))
    db._set_lineage_assignment(hashval2, set([ lin ]))

    assignments = lca_utils.gather_assignments([hashval, hashval2], [db])
    print(assignments)

    assert assignments[hashval] == set([ lin, lin2 ])
    assert assignments[hashval2] == set([ lin ])


def test_count_lca_for_assignments_1():
    # test basic mechanics of gather_assignments function
    hashval = 12345678
    lin = lca_utils.make_lineage('a;b;c')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin ]))

    assignments = lca_utils.gather_assignments([hashval], [db])
    counts = count_lca_for_assignments(assignments)
    print(counts)

    assert len(counts) == 1
    assert counts[lin] == 1


def test_count_lca_for_assignments_2():
    # test basic mechanics of gather_assignments function with two lineages
    hashval = 12345678
    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))

    assignments = lca_utils.gather_assignments([hashval], [db])
    counts = count_lca_for_assignments(assignments)
    print(counts)

    assert counts[lin] == 0
    assert counts[lin2] == 0

    assert len(counts) == 1
    lca_lin = lca_utils.make_lineage('a;b')
    assert counts[lca_lin] == 1


def test_count_lca_for_assignments_3():
    # test basic mechanics of gather_assignments function with two lineages
    # and two hashvals
    hashval = 12345678
    hashval2 = 87654321
    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))
    db._set_lineage_assignment(hashval2, set([ lin ]))

    assignments = lca_utils.gather_assignments([hashval, hashval2], [db])
    counts = count_lca_for_assignments(assignments)
    print(counts)

    assert len(counts) == 2
    assert counts[lin] == 1
    assert counts[lin2] == 0

    lca_lin = lca_utils.make_lineage('a;b')
    assert counts[lca_lin] == 1


def test_count_lca_for_assignments_abund_1():
    # test basic mechanics of gather_assignments function
    hashval = 12345678
    hashval_counts = dict()
    hashval_counts[hashval] = 3

    lin = lca_utils.make_lineage('a;b;c')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin ]))

    assignments = lca_utils.gather_assignments(hashval_counts.keys(), [db])
    counts = count_lca_for_assignments(assignments, hashval_counts)
    print(counts)

    assert len(counts) == 1
    assert counts[lin] == 3


def test_count_lca_for_assignments_abund_2():
    # test basic mechanics of gather_assignments function with two lineages
    hashval = 12345678
    hashval_counts = dict()
    hashval_counts[hashval] = 3

    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))

    assignments = lca_utils.gather_assignments(hashval_counts, [db])
    counts = count_lca_for_assignments(assignments, hashval_counts)
    print(counts)

    assert counts[lin] == 0
    assert counts[lin2] == 0

    assert len(counts) == 1
    lca_lin = lca_utils.make_lineage('a;b')
    assert counts[lca_lin] == 3           # yes!


def test_count_lca_for_assignments_abund_3():
    # test basic mechanics of gather_assignments function with two lineages
    # and two hashvals
    hashval = 12345678
    hashval2 = 87654321
    hashval_counts = dict()
    hashval_counts[hashval] = 2
    hashval_counts[hashval2] = 5

    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))
    db._set_lineage_assignment(hashval2, set([ lin ]))

    assignments = lca_utils.gather_assignments(hashval_counts, [db])
    counts = count_lca_for_assignments(assignments, hashval_counts)
    print(counts)

    assert len(counts) == 2
    assert counts[lin] == 5               # makes sense
    assert counts[lin2] == 0              # makes sense

    lca_lin = lca_utils.make_lineage('a;b')
    assert counts[lca_lin] == 2           # yes!

def test_count_lca_for_assignments_abund_4():
    # test basic mechanics of gather_assignments function with three lineages
    # and three hashvals
    hashval = 12345678
    hashval2 = 87654321
    hashval3 = 34567891
    hashval_counts = dict()
    hashval_counts[hashval] = 2
    hashval_counts[hashval2] = 5
    hashval_counts[hashval3] = 3

    lin = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')
    lin3 = lca_utils.make_lineage('a;b;d;e')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ])) # lca: a;b
    db._set_lineage_assignment(hashval2, set([ lin ])) # lca: a;b;c
    db._set_lineage_assignment(hashval3, set([ lin2, lin3 ])) # a;b;d;e

    assignments = lca_utils.gather_assignments(hashval_counts, [db])
    counts = count_lca_for_assignments(assignments, hashval_counts)
    print(counts)

    assert len(counts) == 3
    assert counts[lin] == 5               # makes sense b/c hashval2
    assert counts[lin2] == 0              # a;b;d (lin2) + a;b;d;e (lin3) -->a;b;d;e (lin3) only
    assert counts[lin3] == 3              # hashval3

    lca_lin = lca_utils.make_lineage('a;b')
    assert counts[lca_lin] == 2           # yes, b/c hashval

def test_count_lca_for_assignments_abund_5():
    # test basic mechanics of gather_assignments function with two lineages
    # and two hashvals when linages match but one has lower taxo detail
    hashval = 12345678
    hashval2 = 87654321
    hashval_counts = dict()
    hashval_counts[hashval] = 2
    hashval_counts[hashval2] = 5

    lin = lca_utils.make_lineage('a;b;d')
    lin2 = lca_utils.make_lineage('a;b;d;e')

    db = FakeLCA_Database()
    db._set_lineage_assignment(hashval, set([ lin, lin2 ]))
    db._set_lineage_assignment(hashval2, set([ lin ]))

    assignments = lca_utils.gather_assignments(hashval_counts, [db])
    counts = count_lca_for_assignments(assignments, hashval_counts)
    print(counts)

    assert len(counts) == 2
    assert counts[lin] == 5               # makes sense
    assert counts[lin2] == 2              # lin+lin2 yield just lin2
