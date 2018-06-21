#! /usr/bin/env python
"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
from __future__ import print_function
import re
import math
import os

from ._minhash import (MinHash, get_minhash_default_seed, get_minhash_max_hash)
DEFAULT_SEED = get_minhash_default_seed()
MAX_HASH = get_minhash_max_hash()

from .due import due, Doi, BibTeX

from .signature import (load_signatures, load_one_signature, SourmashSignature,
                        save_signatures)
from .sbtmh import load_sbt_index, search_sbt_index, create_sbt_index
from . import lca
from . import sbt
from . import sbtmh
from . import sbt_storage
from . import signature

# retrieve VERSION from sourmash/VERSION.
thisdir = os.path.dirname(__file__)
version_file = open(os.path.join(thisdir, 'VERSION'))
VERSION = version_file.read().strip()

due.cite(Doi("10.21105/joss.00027"),
         description="Compute MinHash signatures for nucleotide (DNA/RNA) and protein sequences",
         path="sourmash",
         version=VERSION,
         cite_module=True)

del Doi, BibTeX
