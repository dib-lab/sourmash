from __future__ import print_function
import sys
from setuptools import setup, find_packages
from setuptools import Extension
import os

# retrieve VERSION from sourmash/VERSION.
thisdir = os.path.dirname(__file__)
version_file = open(os.path.join(thisdir, 'sourmash', 'VERSION'))
VERSION = version_file.read().strip()

EXTRA_COMPILE_ARGS = ['-std=c++11', '-pedantic']
EXTRA_LINK_ARGS=[]

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: C++",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

CLASSIFIERS.append("Development Status :: 5 - Production/Stable")

if sys.platform == 'darwin':              # Mac OS X?
    # force 64bit only builds
    EXTRA_COMPILE_ARGS.extend(['-arch', 'x86_64', '-mmacosx-version-min=10.7',
                               '-stdlib=libc++'])

else:                                     # ...likely Linux
   if os.environ.get('SOURMASH_COVERAGE'):
      print('Turning on coverage analysis.')
      EXTRA_COMPILE_ARGS.extend(['-g', '--coverage', '-lgcov'])
      EXTRA_LINK_ARGS.extend(['--coverage', '-lgcov'])
   else:
      EXTRA_COMPILE_ARGS.append('-O3')

SETUP_METADATA = \
               {
    "name": "sourmash",
    "version": VERSION,
    "description": "tools for comparing DNA sequences with MinHash sketches",
    "url": "https://github.com/dib-lab/sourmash",
    "author": "C. Titus Brown",
    "author_email": "titus@idyll.org",
    "license": "BSD 3-clause",
    "packages": find_packages(),
    "entry_points": {'console_scripts': [
        'sourmash = sourmash.__main__:main'
        ]
    },
    "ext_modules": [Extension("sourmash._minhash",
                               sources=["sourmash/_minhash.pyx",
                                        "third-party/smhasher/MurmurHash3.cc"],
                               depends=["sourmash/kmer_min_hash.hh"],
                               include_dirs=["./sourmash",
                                             "./third-party/smhasher/"],
                               language="c++",
                               extra_compile_args=EXTRA_COMPILE_ARGS,
                               extra_link_args=EXTRA_LINK_ARGS)],
    "install_requires": ["screed>=0.9", "ijson", "khmer>=2.1<3.0", "duecredit>=0.6.3"],
    "setup_requires": ['Cython>=0.25.2', "setuptools>=18.0"],
    "extras_require": {
        'test' : ['pytest', 'pytest-cov', 'numpy', 'matplotlib', 'scipy','recommonmark'],
        'demo' : ['jupyter', 'jupyter_client', 'ipython'],
        'doc' : ['sphinx'],
        },
    "include_package_data": True,
    "package_data": {
        "sourmash": ['*.pxd']
    },
    "classifiers": CLASSIFIERS
    }

setup(**SETUP_METADATA)

