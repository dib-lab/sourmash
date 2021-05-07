# Developer information

## Development environment

You can get the latest development branch with:
```
git clone https://github.com/dib-lab/sourmash.git
```
sourmash runs under Python 3.7 and later.

We recommend using `conda` or `Nix` for setting up an environment for developing
new features, running tests and code quality checks.
Here are some suggestions on how to set them up (note: you only need one =])

### Using mamba (conda alternative)

Follow the [installation instructions](https://github.com/conda-forge/miniforge#install) for
installing `mambaforge` (a conda distribution that uses
[`mamba`](https://github.com/TheSnakePit/mamba)
and the [`conda-forge`](https://conda-forge.org/) channel by default).

Once `mamba` is installed, run
```
mamba create -n sourmash_dev tox-conda rust git compilers pandoc
```
to create an environment called `sourmash_dev` containing the programs needed
for development.

To activate the new environment, run
```
conda activate sourmash_dev
```
and proceed to the ["Running tests and checks"](#running-tests-and-checks) section.

### Using Nix

Follow the [installation instructions](https://nixos.org/manual/nix/stable/#chap-installation)
for setting up Nix in your system (Linux or macOS).

Once Nix is installed, run
```
nix-shell
```
to start an environment ready for [running tests and checks](#running-tests-and-checks).

### General instructions

As long as you have `tox` and a Rust compiler available,
you can skip `mamba` or `Nix`.

For Rust, we suggest using `rustup` to install the Rust environment:
```
curl https://sh.rustup.rs -sSf | sh
```
And for `tox` you can run
```
python -m pip install tox
```

We suggest working on sourmash in a virtualenv; e.g. from within the
cloned repository (and after installing `tox` and Rust), you can do:
```
tox -e dev
. .tox/dev/bin/activate
```

Finally, you can also explicitly install all the Python dependencies for sourmash by running
```
pip install -r requirements.txt
```
(but they are already installed in the virtualenv created with `tox -e dev`).

## Updating your developer environment

To update rust to the latest version, use `rustup update`.

To update your Python dependencies to the latest required for sourmash, you can run `pip install -r requirements.txt`.

## Running tests and checks

We use [`tox`](https://tox.readthedocs.io) for managing dependencies and
running tests and checks during development.
`tox -l` lists available tasks.

You can run tests by invoking `make test` in the sourmash directory;
`tox -e py39` will run the Python tests with Python 3.9,
and `cargo test` will run the Rust tests.

## Adding new changes

We use [`pre-commit`](https://pre-commit.com/) to run automatic checks and fixes
when developing sourmash. You can run it with
```
tox -e fix_lint
```
which prints a "hint" at the end of the run with instructions to set it up to
run automatically every time you run `git commit`.

## Automated tests and code coverage calculation

We use [GitHub Actions][2] for continuous integration.

Code coverage can be viewed interactively at [codecov.io][1].

[1]: https://codecov.io/gh/dib-lab/sourmash/
[2]: https://github.com/dib-lab/sourmash/actions

## Code organization

There are three main components in the sourmash repo:
- Python module (in `src/sourmash/`)
- The command-line interface (in `src/sourmash/cli`)
- The Rust core library  (in `src/core`)

`setup.py` has all the configuration to prepare a Python package containing these three components.
First it compiles the Rust core component into a shared library,
which is wrapped by [CFFI] and exposed to the Python module.

[CFFI]: https://cffi.readthedocs.io/

A short description of the high-level files and dirs in the sourmash repo:
```
.
├── benchmarks/         | Benchmarks for the Python module
├── binder/             | mybinder.org configuration
├── data/               | data used for demos
├── doc/                | the documentation rendered in sourmash.bio
├── include/            | C/C++ header files for using core library
├── src/                |
│   ├── core/           | Code for the core library (Rust)
│   └── sourmash/       | The Python module and CLI code
├── tests/              | Tests for the Python module and CLI
├── utils/              |
├── asv.conf.json       | benchmarking config file (for ASV)
├── Cargo.toml          | Rust definition for a workspace
├── CITATION.cff        | Citation info
├── codemeta.json       | Metadata for software discovery
├── CODE_OF_CONDUCT.rst | Code of conduct
├── CONTRIBUTING.md     | Instruction for contributing to development
├── LICENSE             | License for the repo
├── Makefile            | Entry point for most development tasks
├── MANIFEST.in         | Describes what files to add to the Python package
├── matplotlibrc        | Configuration for matplotlib
├── shell.nix           | Nix configuration for creating a dev environment
├── paper.bib           | References in the JOSS paper
├── paper.md            | JOSS paper content
├── pyproject.toml      | Python project definitions (build system and tooling)
├── README.md           | Info to get started
├── requirements.txt    | Python dependencies for development
├── setup.py            | Entry point for Python package setup
├── setup.cfg           | Python package definitions
└── tox.ini             | Configuration for test automation
```

### The Python module (and CLI)

```
src/sourmash
├── cli/                | Command-line parsing, help messages and overall infrastucture
├── command_compute.py  | compute command implementation
├── command_compute.py  | sketch command implementation
├── commands.py         | implementation for other CLI commands
├── compare.py          | Signature comparison functions
├── _compat.py          | Py2/3 compatibility functions
├── exceptions.py       | Mapping from core library errors to Python exceptions
├── fig.py              | Plotting functions
├── index.py            | Index base class and definitions
├── lca/                | LCA index and utility functions
├── logging.py          | Logging functions (notify, error, set_quiet)
├── __main__.py         | Entry point for the CLI
├── _minhash.py         | MinHash sketch implementation (calls the core library)
├── np_utils.py         | NumPy utils
├── sbt*.py             | SBT implementation
├── search.py           | search functions for indices (search, gather)
├── sig                 | signature manipulation functions
│   └── __main__.py     | implementation for `sourmash sig` commands
├── signature_json.py   | signature parsing code (to/from JSON)
├── signature.py        | signature class and methods
├── sourmash_args.py    | convenient shortcuts for CLI usage
└── utils.py            | Convenience functions to interact with core library
```

### The Rust core library

This is completely defined in `src/core` to avoid mixing with the code of other components
(and trying to make it easier to reason about changes).
If you're only working on the core,
you don't need to change any files outside this directory.

This is also published to [crates.io] (the Rust package repository) and [NPM],
after it is compiled to Webassembly.
The GitHub Actions workflow publishes new versions automatically to these
repositories.

[crates.io]: https://crates.io/crates/sourmash
[NPM]: https://www.npmjs.com/package/sourmash

```
src/core
├── benches/             | Benchmarks for the core library
├── Cargo.toml           | Crate definition and metadata
├── cbindgen.toml        | Configuration for cbindgen (the C header generator)
├── examples/            | Examples using the crate API
├── README.md            | Containing links to CI, docs and general info about crate.
├── src                  |
│   ├── cmd.rs           | High-level commands (search, index, compute...)
│   ├── errors.rs        | All the errors generated by this crate
│   ├── ffi/             | FFI-related functions. They are exported to a C header by cbindgen.
│   ├── from.rs          | Conversion methods for other crates
│   ├── index/           | Index methods. An index is a collection of signatures, optimized for searching.
│   ├── lib.rs           | Entry point for the library, control the exposed public API.
│   ├── signature.rs     | Signature methods. A signature is a collection of sketches.
│   ├── sketch/          | Sketch methods. A sketch is compressed representation of data.
│   └── wasm.rs          | Webassembly API.
└── tests/               | Integration tests (using the public API of the crate)
```

### Exposing new functions on the FFI

If you change anything in `src/core/src/ffi` (where the boundary between Rust
and C is defined) you need to regenerate the `include/sourmash.h` header,
and potentially fix any differences in the Python CFFI layer (which reads the C
header file and expose functionality to Python).

To regenerate the C header, run
```
$ make include/sourmash.h
```
This requires a nightly Rust compiler and `cbindgen`.
They can be installed by running
```
$ rustup toolchain add nightly
$ cargo install --force cbindgen
```

### Changing code touching all layers: an example PR

Luiz wrote a [blog post] describing a PR that changes code at the Python API down to the Rust code library,
including some tools for evaluating performance changes.

[blog post]: https://blog.luizirber.org/2020/01/10/sourmash-pr/

## Versioning

We use [`setuptools_scm`] to generate versions based on git tags.
Versions are tagged in a `vMAJOR.MINOR.PATH` format,
following the [Semantic Versioning] convention.
From their definition:

> Given a version number MAJOR.MINOR.PATCH, increment the:
>
> MAJOR version when you make incompatible API changes,
> MINOR version when you add functionality in a backwards compatible manner, and
> PATCH version when you make backwards compatible bug fixes.

[`setuptools_scm`]: https://github.com/pypa/setuptools_scm
[Semantic Versioning]: https://semver.org/

For the Rust core library we use `rMAJOR.MINOR.PATCH`
(note it starts with `r`, and not `v`).
The Rust version is not automated,
and must be bumped in `src/core/Cargo.toml`.

## Common errors and solutions

### Cannot import name `to_bytes` from `sourmash.minhash`

If you are getting an error that contains `ImportError: cannot import name 'to_bytes' from 'sourmash.minhash'`,
then it's likely you need to update Rust and clean up your environment.
Some installation issues can be solved by simply removing the intermediate build files with:

```
make clean
```

## Notes on implementing a new `Index`

- required method: `signatures`
- after that (if there are more efficient impls):
  * `find`
  * `search`/`gather`

## Contents

```{toctree}
:maxdepth: 2

release
requirements
storage
release-notes/releases
```
