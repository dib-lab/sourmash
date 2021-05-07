use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

use crate::signature::SigsTrait;
use crate::sketch::hyperloglog::HyperLogLog;
use crate::traits::Update;

use crate::ffi::minhash::SourmashKmerMinHash;
use crate::ffi::utils::ForeignObject;

pub struct SourmashHyperLogLog;

impl ForeignObject for SourmashHyperLogLog {
    type RustObject = HyperLogLog;
}

#[no_mangle]
pub unsafe extern "C" fn hll_new() -> *mut SourmashHyperLogLog {
    SourmashHyperLogLog::from_rust(HyperLogLog::default())
}

#[no_mangle]
pub unsafe extern "C" fn hll_free(ptr: *mut SourmashHyperLogLog) {
    SourmashHyperLogLog::drop(ptr);
}

ffi_fn! {
unsafe fn hll_with_error_rate(
    error_rate: f64,
    ksize: usize,
) -> Result<*mut SourmashHyperLogLog> {
    let hll = HyperLogLog::with_error_rate(error_rate, ksize)?;
    Ok(SourmashHyperLogLog::from_rust(hll))
}
}

#[no_mangle]
pub unsafe extern "C" fn hll_ksize(ptr: *const SourmashHyperLogLog) -> usize {
    SourmashHyperLogLog::as_rust(ptr).ksize()
}

#[no_mangle]
pub unsafe extern "C" fn hll_cardinality(ptr: *const SourmashHyperLogLog) -> usize {
    SourmashHyperLogLog::as_rust(ptr).cardinality()
}

#[no_mangle]
pub unsafe extern "C" fn hll_similarity(
    ptr: *const SourmashHyperLogLog,
    optr: *const SourmashHyperLogLog,
) -> f64 {
    SourmashHyperLogLog::as_rust(ptr).similarity(SourmashHyperLogLog::as_rust(optr))
}

#[no_mangle]
pub unsafe extern "C" fn hll_containment(
    ptr: *const SourmashHyperLogLog,
    optr: *const SourmashHyperLogLog,
) -> f64 {
    SourmashHyperLogLog::as_rust(ptr).containment(SourmashHyperLogLog::as_rust(optr))
}

#[no_mangle]
pub unsafe extern "C" fn hll_intersection_size(
    ptr: *const SourmashHyperLogLog,
    optr: *const SourmashHyperLogLog,
) -> usize {
    SourmashHyperLogLog::as_rust(ptr).intersection(SourmashHyperLogLog::as_rust(optr))
}

ffi_fn! {
unsafe fn hll_add_sequence(
  ptr: *mut SourmashHyperLogLog,
  sequence: *const c_char,
  insize: usize,
  force: bool
) -> Result<()> {

    let hll = SourmashHyperLogLog::as_rust_mut(ptr);

    let buf = {
        assert!(!ptr.is_null());
        slice::from_raw_parts(sequence as *mut u8, insize)
    };

    hll.add_sequence(buf, force)
}
}

#[no_mangle]
pub unsafe extern "C" fn hll_add_hash(ptr: *mut SourmashHyperLogLog, hash: u64) {
    let hll = SourmashHyperLogLog::as_rust_mut(ptr);
    hll.add_hash(hash);
}

ffi_fn! {
unsafe fn hll_merge(
    ptr: *mut SourmashHyperLogLog,
    optr: *const SourmashHyperLogLog,
) {
    let hll = SourmashHyperLogLog::as_rust_mut(ptr);
    let ohll = SourmashHyperLogLog::as_rust(optr);

    // FIXME raise an exception properly
    hll.merge(ohll)?;
}
}

ffi_fn! {
unsafe fn hll_update_mh(
    ptr: *mut SourmashHyperLogLog,
    optr: *const SourmashKmerMinHash,
) {
    let hll = SourmashHyperLogLog::as_rust_mut(ptr);
    let mh = SourmashKmerMinHash::as_rust(optr);

    mh.update(hll)?
}
}

#[no_mangle]
pub unsafe extern "C" fn hll_matches(
    ptr: *const SourmashHyperLogLog,
    mh_ptr: *const SourmashKmerMinHash,
) -> usize {
    let hll = SourmashHyperLogLog::as_rust(ptr);
    let mh_hll = SourmashKmerMinHash::as_rust(mh_ptr).as_hll();

    hll.intersection(&mh_hll)
}

ffi_fn! {
unsafe fn hll_from_path(filename: *const c_char) -> Result<*mut SourmashHyperLogLog> {
    // FIXME use buffer + len instead of c_str
    let c_str = {
        assert!(!filename.is_null());

        CStr::from_ptr(filename)
    };

    let (mut input, _) = niffler::from_path(c_str.to_str()?)?;
    let hll = HyperLogLog::from_reader(&mut input)?;

    Ok(SourmashHyperLogLog::from_rust(hll))
}
}

ffi_fn! {
unsafe fn hll_from_buffer(ptr: *const c_char, insize: usize) -> Result<*mut SourmashHyperLogLog> {
    // FIXME use SourmashSlice_u8?
    let buf = {
        assert!(!ptr.is_null());
        slice::from_raw_parts(ptr as *mut u8, insize)
    };

    let hll = HyperLogLog::from_reader(buf)?;

    Ok(SourmashHyperLogLog::from_rust(hll))
}
}

ffi_fn! {
unsafe fn hll_save(ptr: *const SourmashHyperLogLog, filename: *const c_char) -> Result<()> {
    let hll = SourmashHyperLogLog::as_rust(ptr);

    // FIXME use buffer + len instead of c_str
    let c_str = {
        assert!(!filename.is_null());

        CStr::from_ptr(filename)
    };

    hll.save(c_str.to_str()?)?;

    Ok(())
}
}

ffi_fn! {
unsafe fn hll_to_buffer(ptr: *const SourmashHyperLogLog, size: *mut usize) -> Result<*const u8> {
    let hll = SourmashHyperLogLog::as_rust(ptr);

    // TODO: remove this
    let compression = 1;

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
      hll.save_to_writer(&mut writer)?;
    }

    let b = buffer.into_boxed_slice();
    *size = b.len();

    // FIXME use SourmashSlice_u8?
    Ok(Box::into_raw(b) as *const u8)
}
}
