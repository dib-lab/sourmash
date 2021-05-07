# -*- coding: UTF-8 -*-

from struct import pack, unpack
import sys
from tempfile import NamedTemporaryFile

from ._lowlevel import ffi, lib
from .minhash import to_bytes, MinHash
from .utils import RustObject, rustcall, decode_str
from .exceptions import SourmashError


class Nodegraph(RustObject):
    __dealloc_func__ = lib.nodegraph_free

    def __init__(self, ksize, starting_size, n_tables):
        self._objptr = lib.nodegraph_with_tables(ksize, int(starting_size), n_tables)

    @staticmethod
    def load(filename):
        ng_ptr = rustcall(lib.nodegraph_from_path, to_bytes(filename))
        return Nodegraph._from_objptr(ng_ptr)

    @staticmethod
    def from_buffer(buf):
        ng_ptr = rustcall(lib.nodegraph_from_buffer, buf, len(buf))
        return Nodegraph._from_objptr(ng_ptr)

    def save(self, filename, version=5):
        assert version >= 4
        if version == 4:
            self._methodcall(lib.nodegraph_save_khmer, to_bytes(filename))
        else:
            self._methodcall(lib.nodegraph_save, to_bytes(filename))

    def to_bytes(self, compression=1):
        size = ffi.new("uintptr_t *")
        rawbuf = self._methodcall(lib.nodegraph_to_buffer, compression, size)
        size = size[0]

        rawbuf = ffi.gc(rawbuf, lambda o: lib.nodegraph_buffer_free(o, size), size)
        buf = ffi.buffer(rawbuf, size)

        return buf

    def update(self, other):
        if isinstance(other, Nodegraph):
            return self._methodcall(lib.nodegraph_update, other._objptr)
        elif isinstance(other, MinHash):
            return self._methodcall(lib.nodegraph_update_mh, other._objptr)
        else:
            # FIXME: we could take sets here too (or anything that can be
            # converted to a list of ints...)
            raise TypeError("Must be a Nodegraph or MinHash")

    def count(self, h):
        if isinstance(h, str):
            return self._methodcall(lib.nodegraph_count_kmer, to_bytes(h))
        return self._methodcall(lib.nodegraph_count, h)

    def get(self, h):
        if isinstance(h, str):
            return self._methodcall(lib.nodegraph_get_kmer, to_bytes(h))
        return self._methodcall(lib.nodegraph_get, h)

    def n_occupied(self):
        return self._methodcall(lib.nodegraph_noccupied)

    def ksize(self):
        return self._methodcall(lib.nodegraph_ksize)

    def hashsizes(self):
        size = ffi.new("uintptr_t *")
        ptr = self._methodcall(lib.nodegraph_hashsizes, size)
        size = size[0]
        hashsizes = ffi.unpack(ptr, size)
        lib.kmerminhash_slice_free(ptr, size)

        return hashsizes

    @property
    def expected_collisions(self):
        return self._methodcall(lib.nodegraph_expected_collisions)

    def matches(self, mh):
        if not isinstance(mh, MinHash):
            # FIXME: we could take sets here too (or anything that can be
            # converted to a list of ints...)
            raise ValueError("mh must be a MinHash")

        return self._methodcall(lib.nodegraph_matches, mh._objptr)

    def to_khmer_nodegraph(self):
        import khmer
        try:
            load_nodegraph = khmer.load_nodegraph
        except AttributeError:
            load_nodegraph = khmer.Nodegraph.load

        with NamedTemporaryFile() as f:
            self.save(f.name, version=4)
            f.file.flush()
            f.file.seek(0)
            return load_nodegraph(f.name)


def extract_nodegraph_info(filename):
    """Open the given nodegraph file and return a tuple of information.

    Returns: the k-mer size, the table size, the number of tables, the version
    of the table format, and the type of table flag.

    Keyword argument:
    filename -- the name of the nodegraph file to inspect
    """
    ksize = None
    n_tables = None
    table_size = None
    signature = None
    version = None
    ht_type = None
    occupied = None

    uint_size = len(pack('I', 0))
    uchar_size = len(pack('B', 0))
    ulonglong_size = len(pack('Q', 0))

    try:
        with open(filename, 'rb') as nodegraph:
            signature, = unpack('4s', nodegraph.read(4))
            version, = unpack('B', nodegraph.read(1))
            ht_type, = unpack('B', nodegraph.read(1))
            ksize, = unpack('I', nodegraph.read(uint_size))
            n_tables, = unpack('B', nodegraph.read(uchar_size))
            occupied, = unpack('Q', nodegraph.read(ulonglong_size))
            table_size, = unpack('Q', nodegraph.read(ulonglong_size))
        if signature != b"OXLI":
            raise ValueError("Node graph '{}' is missing file type "
                             "signature".format(filename) + str(signature))
    except:
        raise ValueError("Node graph '{}' is corrupt ".format(filename))

    return ksize, round(table_size, -2), n_tables, version, ht_type, occupied


def calc_expected_collisions(graph, force=False, max_false_pos=.2):
    fp_all = graph.expected_collisions

    if fp_all > max_false_pos:
        print("**", file=sys.stderr)
        print("** ERROR: the graph structure is too small for ",
              file=sys.stderr)
        print("** this data set.  Increase data structure size.",
              file=sys.stderr)
        print("** Do not use these results!!", file=sys.stderr)
        print("**", file=sys.stderr)
        print("** (estimated false positive rate of %.3f;" % fp_all,
              file=sys.stderr, end=' ')
        print("max recommended %.3f)" % max_false_pos, file=sys.stderr)
        print("**", file=sys.stderr)

        if not force:
            raise SystemExit(1)

    return fp_all
