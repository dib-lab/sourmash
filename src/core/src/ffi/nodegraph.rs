use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

use crate::sketch::nodegraph::Nodegraph;
use crate::traits::Update;

use crate::ffi::minhash::SourmashKmerMinHash;
use crate::ffi::utils::ForeignObject;

pub struct SourmashNodegraph;

impl ForeignObject for SourmashNodegraph {
    type RustObject = Nodegraph;
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_new() -> *mut SourmashNodegraph {
    SourmashNodegraph::from_rust(Nodegraph::default())
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_free(ptr: *mut SourmashNodegraph) {
    SourmashNodegraph::drop(ptr);
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_buffer_free(ptr: *mut u8, insize: usize) {
    if ptr.is_null() {
        return;
    }
    Vec::from_raw_parts(ptr as *mut u8, insize, insize);
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_with_tables(
    ksize: usize,
    starting_size: usize,
    n_tables: usize,
) -> *mut SourmashNodegraph {
    let ng = Nodegraph::with_tables(starting_size, n_tables, ksize);
    SourmashNodegraph::from_rust(ng)
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_count(ptr: *mut SourmashNodegraph, h: u64) -> bool {
    let ng = SourmashNodegraph::as_rust_mut(ptr);
    ng.count(h)
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_count_kmer(
    ptr: *mut SourmashNodegraph,
    kmer: *const c_char,
) -> bool {
    let ng = SourmashNodegraph::as_rust_mut(ptr);

    // FIXME use buffer + len instead of cstr
    let c_str = {
        assert!(!kmer.is_null());

        CStr::from_ptr(kmer)
    };

    ng.count_kmer(c_str.to_bytes())
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_get(ptr: *const SourmashNodegraph, h: u64) -> usize {
    let ng = SourmashNodegraph::as_rust(ptr);
    ng.get(h)
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_get_kmer(
    ptr: *const SourmashNodegraph,
    kmer: *const c_char,
) -> usize {
    let ng = SourmashNodegraph::as_rust(ptr);

    // FIXME use buffer + len instead of cstr
    let c_str = {
        assert!(!kmer.is_null());

        CStr::from_ptr(kmer)
    };

    ng.get_kmer(c_str.to_bytes())
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_expected_collisions(ptr: *const SourmashNodegraph) -> f64 {
    let ng = SourmashNodegraph::as_rust(ptr);
    ng.expected_collisions()
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_ksize(ptr: *const SourmashNodegraph) -> usize {
    let ng = SourmashNodegraph::as_rust(ptr);
    ng.ksize()
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_hashsizes(
    ptr: *const SourmashNodegraph,
    size: *mut usize,
) -> *const u64 {
    let ng = SourmashNodegraph::as_rust(ptr);
    let st = ng.tablesizes();

    let b = st.into_boxed_slice();
    *size = b.len();

    // FIXME: Use SourmashSlice_u64?
    Box::into_raw(b) as *const u64
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_ntables(ptr: *const SourmashNodegraph) -> usize {
    let ng = SourmashNodegraph::as_rust(ptr);
    ng.ntables()
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_noccupied(ptr: *const SourmashNodegraph) -> usize {
    let ng = SourmashNodegraph::as_rust(ptr);
    ng.noccupied()
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_matches(
    ptr: *const SourmashNodegraph,
    mh_ptr: *const SourmashKmerMinHash,
) -> usize {
    let ng = SourmashNodegraph::as_rust(ptr);
    let mh = SourmashKmerMinHash::as_rust(mh_ptr);
    ng.matches(mh)
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_update(
    ptr: *mut SourmashNodegraph,
    optr: *const SourmashNodegraph,
) {
    let ng = SourmashNodegraph::as_rust_mut(ptr);
    let ong = SourmashNodegraph::as_rust(optr);

    // FIXME raise an exception properly
    ong.update(ng).unwrap();
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_update_mh(
    ptr: *mut SourmashNodegraph,
    optr: *const SourmashKmerMinHash,
) {
    let ng = SourmashNodegraph::as_rust_mut(ptr);
    let mh = SourmashKmerMinHash::as_rust(optr);

    mh.update(ng).unwrap();
}

ffi_fn! {
unsafe fn nodegraph_from_path(filename: *const c_char) -> Result<*mut SourmashNodegraph> {
    // FIXME use buffer + len instead of c_str
    let c_str = {
        assert!(!filename.is_null());

        CStr::from_ptr(filename)
    };

    let (mut input, _) = niffler::from_path(c_str.to_str()?)?;
    let ng = Nodegraph::from_reader(&mut input)?;

    Ok(SourmashNodegraph::from_rust(ng))
}
}

ffi_fn! {
unsafe fn nodegraph_from_buffer(ptr: *const c_char, insize: usize) -> Result<*mut SourmashNodegraph> {
    // FIXME use SourmashSlice_u8?
    let buf = {
        assert!(!ptr.is_null());
        slice::from_raw_parts(ptr as *mut u8, insize)
    };

    let ng = Nodegraph::from_reader(buf)?;

    Ok(SourmashNodegraph::from_rust(ng))
}
}

ffi_fn! {
unsafe fn nodegraph_save(ptr: *const SourmashNodegraph, filename: *const c_char) -> Result<()> {
    let ng = SourmashNodegraph::as_rust(ptr);

    // FIXME use buffer + len instead of c_str
    let c_str = {
        assert!(!filename.is_null());

        CStr::from_ptr(filename)
    };

    ng.save(c_str.to_str()?)?;

    Ok(())
}
}

ffi_fn! {
unsafe fn nodegraph_to_buffer(ptr: *const SourmashNodegraph, compression: u8, size: *mut usize) -> Result<*const u8> {
    let ng = SourmashNodegraph::as_rust(ptr);

    let mut buffer = vec![];
    {
      let mut writer = if compression > 0 {
          let level = match compression {
            1 => niffler::compression::Level::One,
            2 => niffler::compression::Level::Two,
            3 => niffler::compression::Level::Three,
            4 => niffler::compression::Level::Four,
            5 => niffler::compression::Level::Five,
            6 => niffler::compression::Level::Six,
            7 => niffler::compression::Level::Seven,
            8 => niffler::compression::Level::Eight,
            _ => niffler::compression::Level::Nine,
          };

          niffler::get_writer(Box::new(&mut buffer),
                              niffler::compression::Format::Gzip,
                              level)?
      } else {
          Box::new(&mut buffer)
      };
      ng.save_to_writer(&mut writer)?;
    }

    let b = buffer.into_boxed_slice();
    *size = b.len();

    // FIXME use SourmashSlice_u8?
    Ok(Box::into_raw(b) as *const u8)
}
}
