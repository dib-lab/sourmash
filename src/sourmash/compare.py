"""Functionality for comparing many signatures, used in sourmash compare."""

import itertools
from functools import partial
import time
import multiprocessing

from .logging import notify
from sourmash.np_utils import to_memmap


def compare_serial(siglist, ignore_abundance, downsample=False):
    """Compare all combinations of signatures and return a matrix
    of similarities. Processes combinations serially on a single
    process. Best to use when there is few signatures.

    :param list siglist: list of signatures to compare
    :param boolean ignore_abundance
        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate the angular
        similarity.
    :param boolean downsample by scaled if True
    :return: np.array similarity matrix
    """
    import numpy as np

    n = len(siglist)

    # Combinations makes all unique sets of pairs, e.g. (A, B) but not (B, A)
    iterator = itertools.combinations(range(n), 2)

    similarities = np.ones((n, n))

    for i, j in iterator:
        similarities[i][j] = similarities[j][i] = siglist[i].similarity(siglist[j], ignore_abundance, downsample)

    return similarities


def compare_serial_containment(siglist, downsample=False):
    """Compare all combinations of signatures and return a matrix
    of containments. Processes combinations serially on a single
    process. Best to only use when there are few signatures.

    :param list siglist: list of signatures to compare
    :param boolean downsample by scaled if True
    :return: np.array similarity matrix
    """
    import numpy as np

    n = len(siglist)

    containments = np.ones((n, n))
    for i in range(n):
        for j in range(n):
            containments[i][j] = siglist[j].contained_by(siglist[i],
                                                         downsample=downsample)

    return containments


def compare_serial_max_containment(siglist, downsample=False):
    """Compare all combinations of signatures and return a matrix
    of max_containments. Processes combinations serially on a single
    process. Best to only use when there are few signatures.

    :param list siglist: list of signatures to compare
    :param boolean downsample by scaled if True
    :return: np.array similarity matrix
    """
    import numpy as np

    n = len(siglist)

    containments = np.ones((n, n))
    for i in range(n):
        for j in range(n):
            containments[i][j] = siglist[j].max_containment(siglist[i],
                                                         downsample=downsample)

    return containments


def similarity_args_unpack(args, ignore_abundance, downsample):
    """Helper function to unpack the arguments. Written to use in pool.imap
    as it can only be given one argument."""
    sig1, sig2 = args
    return sig1.similarity(sig2,
                           ignore_abundance=ignore_abundance,
                           downsample=downsample)


def compare_parallel(siglist, ignore_abundance, downsample, n_jobs):
    """Compare all combinations of signatures and return a matrix
    of similarities. Processes combinations parallely on number of processes
    given by n_jobs

    :param list siglist: list of signatures to compare
    :param boolean ignore_abundance
        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate the angular
        similarity.
    :param boolean downsample by scaled if True
    :param int n_jobs number of processes to run the similarity calculations on
    :return: np.array similarity matrix
    """
    import numpy as np

    # Starting time - calculate time to keep track in case of lengthy siglist
    start_initial = time.time()

    # Create a memory map of the siglist using numpy to avoid memory burden
    # while accessing small parts in it
    siglist, _ = to_memmap(np.array(siglist))
    notify("Created memmapped siglist")

    # Check that length of combinations can result in a square similarity matrix
    length_siglist = len(siglist)

    # Initialize with ones in the diagonal as the similarity of a signature with
    # itself is one
    similarities = np.eye(length_siglist, dtype=np.float64)
    memmap_similarities, filename = to_memmap(similarities)
    notify("Initialized memmapped similarities matrix")

    # To avoid sharing siglist via pickle declaring a function inside function

    def get_similarities_at_index(index):
        """Returns similarities of all the combinations of signature at index in
        the siglist with the rest of the indices starting at index + 1. Doesn't
        redundantly calculate signatures with all the other indices prior to
        index - 1

        :param int index: generate masks from this image
        :param boolean ignore_abundance
            If the sketches are not abundance weighted, or ignore_abundance=True,
            compute Jaccard similarity.

            If the sketches are abundance weighted, calculate the angular
            similarity.
        :param boolean downsample by scaled if True
        :return: list of similarities for the combinations of signature at index
            with rest of the signatures from index+1
        """
        startt = time.time()
        sig_iterator = itertools.product([siglist[index]], siglist[index + 1:])
        func = partial(similarity_args_unpack,
                       ignore_abundance=ignore_abundance,
                       downsample=downsample)
        similarity_list = list(map(func, sig_iterator))
        notify(
            "comparison for index {} done in {:.5f} seconds",
            index,
            time.time() - startt,
            end='\r')
        return similarity_list

    # Initialize multiprocess.pool
    pool = multiprocessing.Pool(processes=n_jobs)

    # Calculate chunk size, by default pool.imap chunk size is 1
    chunksize, extra = divmod(length_siglist, n_jobs)
    if extra:
        chunksize += 1
    notify("Calculated chunk size for multiprocessing")

    # This will not generate the results yet, since pool.imap returns a generator
    result = pool.imap(get_similarities_at_index, range(length_siglist), chunksize=chunksize)
    notify("Initialized multiprocessing pool.imap")

    # Enumerate and calculate similarities at each of the indices
    # and set the results at the appropriate combination coordinate
    # locations inside the similarity matrix
    for index, l in enumerate(result):
        startt = time.time()
        col_idx = index + 1
        for idx_condensed, item in enumerate(l):
            memmap_similarities[index, col_idx + idx_condensed] = memmap_similarities[idx_condensed + col_idx, index] = item
        notify(
            "Setting similarities matrix for index {} done in {:.5f} seconds",
            index,
            time.time() - startt,
            end='\r')
    notify("Setting similarities completed")

    pool.close()
    pool.join()

    notify("Time taken to compare all pairs parallely is {:.5f} seconds ", time.time() - start_initial)
    return np.memmap(filename, dtype=np.float64, shape=(length_siglist, length_siglist))


def compare_all_pairs(siglist, ignore_abundance, downsample=False, n_jobs=None):
    """Compare all combinations of signatures and return a matrix
    of similarities. Processes combinations either serially or
    based on parallely on number of processes given by n_jobs

    :param list siglist: list of signatures to compare
    :param boolean ignore_abundance
        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate the angular
        similarity.
    :param boolean downsample by scaled if True
    :param int n_jobs number of processes to run the similarity calculations on,
    if number of jobs is None or 1, compare serially, otherwise parallely.
    :return: np.array similarity matrix
    """
    if n_jobs is None or n_jobs == 1:
        similarities = compare_serial(siglist, ignore_abundance, downsample)
    else:
        similarities = compare_parallel(siglist, ignore_abundance, downsample, n_jobs)
    return similarities
