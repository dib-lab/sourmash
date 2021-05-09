#! /usr/bin/env python
"""
Save and load MinHash sketches in a JSON format, along with some metadata.
"""
import sys
import os
import weakref
from enum import Enum

from .logging import error
from . import MinHash
from .minhash import to_bytes, FrozenMinHash
from ._lowlevel import ffi, lib
from .utils import RustObject, rustcall, decode_str


SIGNATURE_VERSION = 0.4


class SigInput(Enum):
    FILE_LIKE = 1
    PATH = 2
    BUFFER = 3
    UNKNOWN = 4


class SourmashSignature(RustObject):
    "Main class for signature information."

    __dealloc_func__ = lib.signature_free

    def __init__(self, minhash, name="", filename=""):
        self._objptr = lib.signature_new()

        if name:
            self.name = name
        if filename:
            self.filename = filename

        self.minhash = minhash


    @property
    def minhash(self):
        return FrozenMinHash._from_objptr(
            self._methodcall(lib.signature_first_mh)
        )

    @minhash.setter
    def minhash(self, value):
        # TODO: validate value is a MinHash
        self._methodcall(lib.signature_set_mh, value._objptr)

    def __hash__(self):
        return hash(self.md5sum())

    def __str__(self):
        return self._display_name()

    def __repr__(self):
        name = self.name
        md5pref = self.md5sum()[:8]
        if name == md5pref:
            return "SourmashSignature({})".format(md5pref)
        else: # name != md5pref:
            return "SourmashSignature('{}', {})".format(name, md5pref)

    #def minhashes(self):
    #    size = ffi.new("uintptr_t *")
    #    mhs_ptr = self._methodcall(lib.signature_get_mhs, size)
    #    size = ffi.unpack(size, 1)[0]
    #
    #    mhs = []
    #    for i in range(size):
    #        mh = MinHash._from_objptr(mhs_ptr[i])
    #        mhs.append(mh)
    #
    #    return mhs

    def md5sum(self):
        "Calculate md5 hash of the bottom sketch, specifically."
        return decode_str(self.minhash._methodcall(lib.kmerminhash_md5sum))

    def __eq__(self, other):
        return self._methodcall(lib.signature_eq, other._objptr)

    @property
    def _name(self):
        return decode_str(self._methodcall(lib.signature_get_name))

    @_name.setter
    def _name(self, value):
        self._methodcall(lib.signature_set_name, to_bytes(value))

    @property
    def name(self):
        return decode_str(self._methodcall(lib.signature_get_name))

    @name.setter
    def name(self, value):
        self._methodcall(lib.signature_set_name, to_bytes(value))

    def __ne__(self, other):
        return not self == other

    @property
    def filename(self):
        return decode_str(self._methodcall(lib.signature_get_filename))

    @filename.setter
    def filename(self, value):
        self._methodcall(lib.signature_set_filename, to_bytes(value))

    @property
    def license(self):
        return decode_str(self._methodcall(lib.signature_get_license))

    def _display_name(self, max_length=0):
        name = self._name
        filename = self.filename

        if name:
            if max_length and len(name) > max_length:
                name = name[: max_length - 3] + "..."
        elif filename:
            name = filename
            if max_length and len(name) > max_length:
                name = "..." + name[-max_length + 3 :]
        else:
            name = self.md5sum()[:8]
        assert not max_length or len(name) <= max_length
        return name

    def similarity(self, other, ignore_abundance=False, downsample=False):
        "Compute similarity with the other signature."
        return self.minhash.similarity(other.minhash,
                                       ignore_abundance=ignore_abundance,
                                       downsample=downsample)

    def jaccard(self, other):
        "Compute Jaccard similarity with the other MinHash signature."
        return self.minhash.similarity(other.minhash, ignore_abundance=True,
                                       downsample=False)

    def contained_by(self, other, downsample=False):
        "Compute containment by the other signature. Note: ignores abundance."
        return self.minhash.contained_by(other.minhash, downsample)

    def max_containment(self, other, downsample=False):
        "Compute max containment w/other signature. Note: ignores abundance."
        return self.minhash.max_containment(other.minhash, downsample)

    def add_sequence(self, sequence, force=False):
        self._methodcall(lib.signature_add_sequence, to_bytes(sequence), force)

    def add_protein(self, sequence):
        self._methodcall(lib.signature_add_protein, to_bytes(sequence))

    @staticmethod
    def from_params(params):
        ptr = rustcall(lib.signature_from_params, params._get_objptr())
        return SourmashSignature._from_objptr(ptr)

    def __len__(self):
        return self._methodcall(lib.signature_len)

    def __getstate__(self):  # enable pickling
        return (
            self.minhash,
            self.name,
            self.filename,
        )

    def __setstate__(self, tup):
        (mh, name, filename) = tup
        self.__del__()

        self._objptr = lib.signature_new()
        if name:
            self.name = name
        if filename:
            self.filename = filename
        self.minhash = mh

    def __reduce__(self):
        return (
            SourmashSignature,
            (
                self.minhash,
                self.name,
                self.filename
            ),
        )


def _detect_input_type(data):
    """\
    Determine how to load input from `data`. Returns SigInput enum.

    Checks for:
     - Python file-like objects
     - JSON text (uncompressed sigs)
     - Compressed memory buffers
     - filename
    """
    if hasattr(data, 'read') or hasattr(data, "fileno") or hasattr(data, "mode"):  # file-like object
        return SigInput.FILE_LIKE
    elif hasattr(data, "find"):  # check if it is uncompressed sig
        try:
            if data.find("sourmash_signature") > 0:
                return SigInput.BUFFER
        except TypeError:
            if data.find(b"sourmash_signature") > 0:
                return SigInput.BUFFER
            elif data.startswith(b'\x1F\x8B'):  # gzip compressed
                return SigInput.BUFFER

    try:
        if os.path.exists(data):  # filename
            return SigInput.PATH
    except (ValueError, TypeError):  # No idea...
        return SigInput.UNKNOWN

    return SigInput.UNKNOWN


def load_signatures(
    data, ksize=None, select_moltype=None, ignore_md5sum=False, do_raise=False,
):
    """Load a JSON string with signatures into classes.

    Returns list of SourmashSignature objects.

    Note, the order is not necessarily the same as what is in the source file.
    """
    if ksize is not None:
        ksize = int(ksize)
    else:
        ksize = 0

    if not data:
        return

    if select_moltype is None:
        select_moltype = ffi.NULL
    else:
        try:
            select_moltype = select_moltype.encode("utf-8")
        except AttributeError:
            pass

    input_type = _detect_input_type(data)
    if input_type == SigInput.UNKNOWN:
        if do_raise:
            raise ValueError("Error in parsing signature; quitting. Cannot open file or invalid signature")
        return

    size = ffi.new("uintptr_t *")

    try:
        if input_type == SigInput.FILE_LIKE:
            if hasattr(data, "mode") and "t" in data.mode:  # need to reopen handler as binary
                data = data.buffer

            buf = data.read()
            data.close()
            data = buf
            input_type = SigInput.BUFFER

        elif input_type == SigInput.PATH:
            sigs_ptr = rustcall(
                lib.signatures_load_path,
                data.encode("utf-8"),
                ignore_md5sum,
                ksize,
                select_moltype,
                size,
            )

        if input_type == SigInput.BUFFER:
            if hasattr(data, "encode"):
                data = data.encode("utf-8")

            sigs_ptr = rustcall(
                lib.signatures_load_buffer,
                data,
                len(data),
                ignore_md5sum,
                ksize,
                select_moltype,
                size,
            )

        size = size[0]

        sigs = []
        for i in range(size):
            sig = SourmashSignature._from_objptr(sigs_ptr[i])
            sigs.append(sig)

        for sig in sigs:
            yield sig

    except Exception as e:
        if do_raise:
            raise


def load_one_signature(data, ksize=None, select_moltype=None, ignore_md5sum=False):
    sigiter = load_signatures(
        data, ksize=ksize, select_moltype=select_moltype, ignore_md5sum=ignore_md5sum
    )

    try:
        first_sig = next(sigiter)
    except StopIteration:
        raise ValueError("no signatures to load")

    try:
        next(sigiter)
    except StopIteration:
        return first_sig

    raise ValueError("expected to load exactly one signature")


def save_signatures(siglist, fp=None, compression=0):
    "Save multiple signatures into a JSON string (or into file handle 'fp')"
    attached_refs = weakref.WeakKeyDictionary()

    # get list of rust objects
    collected = []
    for obj in siglist:
        rv = obj._get_objptr()
        attached_refs[rv] = obj
        collected.append(rv)
    siglist_c = ffi.new("SourmashSignature*[]", collected)

    size = ffi.new("uintptr_t *")

    # save signature into a string (potentially compressed)
    rawbuf = rustcall(lib.signatures_save_buffer, siglist_c, len(collected),
                      compression, size)
    size = size[0]

    # associate a finalizer with rawbuf so that it gets freed
    buf = ffi.gc(rawbuf, lambda o: lib.nodegraph_buffer_free(o, size), size)
    if compression:
        result = ffi.buffer(buf, size)[:]
    else:
        result = ffi.string(buf, size)

    if fp is None:                        # return string
        return result
    else:
        try:                              # write to file
            fp.write(result)
        except TypeError:
            fp.write(result.decode('utf-8'))
        return None
