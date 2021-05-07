import os
import sys

from setuptools import setup


DEBUG_BUILD = os.environ.get("SOURMASH_DEBUG") == "1"


def build_native(spec):
    cmd = ["cargo", "build",
           "--manifest-path", "src/core/Cargo.toml",
           "--features", "experimental,parallel",
           "--lib"]

    target = "debug"
    if not DEBUG_BUILD:
        cmd.append("--release")
        target = "release"

    build = spec.add_external_build(cmd=cmd, path=".")

    rtld_flags = ["NOW"]
    if sys.platform == "darwin":
        rtld_flags.append("NODELETE")
    spec.add_cffi_module(
        module_path="sourmash._lowlevel",
        dylib=lambda: build.find_dylib("sourmash", in_path="target/%s" % target),
        header_filename=lambda: build.find_header("sourmash.h", in_path="include"),
        rtld_flags=rtld_flags,
    )

setup(
  milksnake_tasks=[build_native],
  package_dir={"": "src"},
)
