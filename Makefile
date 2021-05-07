PYTHON ?= python

all: build

.PHONY:

build: .PHONY
	$(PYTHON) setup.py build_ext -i

clean:
	$(PYTHON) setup.py clean --all
	rm -f src/sourmash/*.so
	cd doc && make clean

install: all
	$(PYTHON) setup.py install

dist: FORCE
	$(PYTHON) setup.py sdist

test:
	tox -e py38
	cargo test

doc: .PHONY
	tox -e docs

include/sourmash.h: src/core/src/lib.rs \
                    src/core/src/ffi/hyperloglog.rs \
                    src/core/src/ffi/minhash.rs \
                    src/core/src/ffi/signature.rs \
                    src/core/src/ffi/nodegraph.rs \
                    src/core/src/ffi/index/mod.rs \
                    src/core/src/ffi/index/revindex.rs \
                    src/core/src/errors.rs
	cd src/core && \
	RUSTUP_TOOLCHAIN=nightly cbindgen -c cbindgen.toml . -o ../../$@

coverage: all
	tox -e coverage

benchmark:
	tox -e asv
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
