from __future__ import print_function, unicode_literals

import os
import shutil

import pytest

from sourmash_lib import signature
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import SigLeaf, search_minhashes
from sourmash_lib.sbt_storage import IPFSStorage, FSStorage

from . import sourmash_tst_utils as utils


def test_sbt_ipfsstorage(tmpdir):
    pytest.importorskip("ipfsapi")

    factory = GraphFactory(31, 1e5, 4)

    tree = SBT(factory)

    for f in utils.SIG_FILES:
        sig = next(signature.load_signatures(utils.get_test_data(f)))
        leaf = SigLeaf(os.path.basename(f), sig)
        tree.add_node(leaf)
        to_search = leaf

    print("*" * 60)
    print("{}:".format(to_search.metadata))
    old_result = {str(s) for s in tree.find(search_minhashes, to_search.data, 0.1)}
    print(*old_result, sep="\n")

    try:
        with IPFSStorage() as storage:
            tree.save(str(tmpdir.join("tree")), storage=storage)
    except NotImplementedError:
        pytest.xfail("Using a Read-only client for IPFS")

    with IPFSStorage() as storage:
        tree = SBT.load(
            str(tmpdir.join("tree")), leaf_loader=SigLeaf.load, storage=storage
        )

        print("*" * 60)
        print("{}:".format(to_search.metadata))
        new_result = {str(s) for s in tree.find(search_minhashes, to_search.data, 0.1)}
        print(*new_result, sep="\n")

        assert old_result == new_result


def test_storage_convert(tmpdir):
    testdata = utils.get_test_data("v2.sbt.json")
    testsbt = tmpdir.join("v2.sbt.json")

    shutil.copyfile(testdata, str(testsbt))
    shutil.copytree(
        os.path.join(os.path.dirname(testdata), ".sbt.v2"), str(tmpdir.join(".sbt.v2"))
    )

    original = SBT.load(str(testsbt), leaf_loader=SigLeaf.load)

    args = ["storage", "convert", "-b", "ipfs", testsbt]
    status, out, err = utils.runscript(
        "sourmash", args, in_directory=str(tmpdir), fail_ok=True
    )
    if not status and "ipfs.exceptions.ConnectionError" in err:
        raise pytest.xfail("ipfs probably not running")

    ipfs = SBT.load(str(testsbt), leaf_loader=SigLeaf.load)

    assert len(original) == len(ipfs)
    assert all(n1[1].name == n2[1].name for (n1, n2) in zip(original, ipfs))

    args = [
        "storage",
        "convert",
        "-b",
        """'TarStorage("{}")'""".format(tmpdir.join("v2.sbt.tar.gz")),
        str(testsbt),
    ]
    status, out, err = utils.runscript("sourmash", args, in_directory=str(tmpdir))
    tar = SBT.load(str(testsbt), leaf_loader=SigLeaf.load)

    assert len(original) == len(tar)
    assert all(n1[1].name == n2[1].name for (n1, n2) in zip(original, tar))


def test_prepare_index(tmpdir):
    pytest.importorskip("ipfsapi")

    try:
        with IPFSStorage() as storage:
            for f in utils.SIG_FILES:
                with open(utils.get_test_data(f), 'rb') as data:
                    storage.save('', data.read())
    except NotImplementedError:
        pytest.xfail("Using a Read-only client for IPFS")

    testdata = utils.get_test_data("ipfs_leaves.sbt.json")
    testsbt = tmpdir.join("tree.sbt.json")

    shutil.copyfile(testdata, str(testsbt))

    args = [
        "prepare",
        str(testsbt),
    ]
    status, out, err = utils.runscript("sourmash", args, in_directory=str(tmpdir))
    prepared_sbt = SBT.load(str(testsbt), leaf_loader=SigLeaf.load)
    assert not(isinstance(prepared_sbt.storage, IPFSStorage))
    assert isinstance(prepared_sbt.storage, FSStorage)

    sig = utils.get_test_data(utils.SIG_FILES[0])
    status, out, err = utils.runscript('sourmash',
                                       ['search', sig, str(testsbt)],
                                       in_directory=str(tmpdir))
    assert '3 matches:' in out
