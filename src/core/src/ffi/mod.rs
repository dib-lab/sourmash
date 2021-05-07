//! # Foreign Function Interface for calling sourmash from a C API
//!
//! Primary client for now is the Python version, using CFFI and milksnake.
#![allow(clippy::missing_safety_doc)]

#[macro_use]
pub mod utils;

pub mod cmd;
pub mod hyperloglog;
pub mod index;
pub mod minhash;
pub mod nodegraph;
pub mod signature;

use std::ffi::CStr;
use std::os::raw::c_char;

use crate::_hash_murmur;

#[no_mangle]
pub unsafe extern "C" fn hash_murmur(kmer: *const c_char, seed: u64) -> u64 {
    let c_str = {
        assert!(!kmer.is_null());

        CStr::from_ptr(kmer)
    };

    _hash_murmur(c_str.to_bytes(), seed)
}
