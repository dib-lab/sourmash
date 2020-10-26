PYTHON ?= python

all: build

.PHONY:

build: .PHONY
	$(PYTHON) setup.py build_ext -i

clean:
	$(PYTHON) setup.py clean --all
	rm -f sourmash/*.so
	cd doc && make clean

install: all
	$(PYTHON) setup.py install

dist: FORCE
	$(PYTHON) setup.py sdist

test: all
	$(PYTHON) -m pip install -e '.[test]'
	$(PYTHON) -m pytest
	cargo test

doc: build .PHONY
	cd doc && make html

include/sourmash.h: src/core/src/lib.rs \
                    src/core/src/ffi/minhash.rs \
                    src/core/src/ffi/signature.rs \
                    src/core/src/ffi/nodegraph.rs \
                    src/core/src/ffi/lca_db.rs \
                    src/core/src/errors.rs
	cd src/core && \
	RUSTUP_TOOLCHAIN=nightly cbindgen -c cbindgen.toml . -o ../../$@

coverage: all
	$(PYTHON) setup.py build_ext -i
	$(PYTHON) -m pytest --cov=. --cov-report term-missing

benchmark:
	asv continuous latest `git rev-parse HEAD`
	cargo bench

check:
	cargo build
	cargo test
	cargo bench

last-tag:
	git fetch -p -q; git tag -l | sort -V | tail -1

wasm:
	wasm-pack build src/core -d ../../pkg

wasi:
	cargo wasi build

FORCE:
