use std::convert::TryInto;
use std::ffi::CStr;
use std::io;
use std::os::raw::c_char;
use std::slice;

use crate::encodings::HashFunctions;
use crate::signature::Signature;
use crate::sketch::Sketch;

use crate::ffi::cmd::compute::SourmashComputeParameters;
use crate::ffi::minhash::SourmashKmerMinHash;
use crate::ffi::utils::{ForeignObject, SourmashStr};

pub struct SourmashSignature;

impl ForeignObject for SourmashSignature {
    type RustObject = Signature;
}

// Signature methods

#[no_mangle]
pub unsafe extern "C" fn signature_new() -> *mut SourmashSignature {
    SourmashSignature::from_rust(Signature::default())
}

#[no_mangle]
pub unsafe extern "C" fn signature_from_params(
    ptr: *const SourmashComputeParameters,
) -> *mut SourmashSignature {
    let params = SourmashComputeParameters::as_rust(ptr);
    SourmashSignature::from_rust(Signature::from_params(params))
}

#[no_mangle]
pub unsafe extern "C" fn signature_free(ptr: *mut SourmashSignature) {
    SourmashSignature::drop(ptr);
}

#[no_mangle]
pub unsafe extern "C" fn signature_len(ptr: *const SourmashSignature) -> usize {
    let sig = SourmashSignature::as_rust(ptr);
    sig.size()
}

ffi_fn! {
unsafe fn signature_add_sequence(ptr: *mut SourmashSignature, sequence: *const c_char, force: bool) ->
    Result<()> {
    let sig = SourmashSignature::as_rust_mut(ptr);

    // FIXME replace with buffer + len
    let c_str = {
        assert!(!sequence.is_null());

        CStr::from_ptr(sequence)
    };

    sig.add_sequence(c_str.to_bytes(), force)
}
}

ffi_fn! {
unsafe fn signature_add_protein(ptr: *mut SourmashSignature, sequence: *const c_char) ->
    Result<()> {
    let sig = SourmashSignature::as_rust_mut(ptr);

    // FIXME replace with buffer + len
    let c_str = {
        assert!(!sequence.is_null());

        CStr::from_ptr(sequence)
    };

    sig.add_protein(c_str.to_bytes())
}
}

ffi_fn! {
unsafe fn signature_set_name(ptr: *mut SourmashSignature, name: *const c_char) ->
    Result<()> {
    let sig = SourmashSignature::as_rust_mut(ptr);

    // FIXME replace with buffer + len
    let c_str = {
        assert!(!name.is_null());

        CStr::from_ptr(name)
    };

    if let Ok(name) = c_str.to_str() {
        sig.set_name(name)
    }
    Ok(())
}
}

ffi_fn! {
unsafe fn signature_set_filename(ptr: *mut SourmashSignature, name: *const c_char) ->
    Result<()> {
    let sig = SourmashSignature::as_rust_mut(ptr);

    // FIXME replace with buffer + len
    let c_str = {
        assert!(!name.is_null());

        CStr::from_ptr(name)
    };

    if let Ok(name) = c_str.to_str() {
        sig.set_filename(name)
    }
    Ok(())
}
}

ffi_fn! {
unsafe fn signature_push_mh(ptr: *mut SourmashSignature, other: *const SourmashKmerMinHash) ->
    Result<()> {
    let sig = SourmashSignature::as_rust_mut(ptr);
    let mh = SourmashKmerMinHash::as_rust(other);
    sig.push(Sketch::MinHash(mh.clone()));
    Ok(())
}
}

ffi_fn! {
unsafe fn signature_set_mh(ptr: *mut SourmashSignature, other: *const SourmashKmerMinHash) ->
    Result<()> {
    let sig = SourmashSignature::as_rust_mut(ptr);
    let mh = SourmashKmerMinHash::as_rust(other);
    sig.reset_sketches();
    sig.push(Sketch::MinHash(mh.clone()));
    Ok(())
}
}

ffi_fn! {
unsafe fn signature_get_name(ptr: *const SourmashSignature) -> Result<SourmashStr> {
    let sig = SourmashSignature::as_rust(ptr);

    if let Some(ref name) = sig.name {
        Ok(name.clone().into())
    } else {
        Ok("".into())
    }
}
}

ffi_fn! {
unsafe fn signature_get_filename(ptr: *const SourmashSignature) -> Result<SourmashStr> {
    let sig = SourmashSignature::as_rust(ptr);
    Ok(sig.filename().into())
}
}

ffi_fn! {
unsafe fn signature_get_license(ptr: *const SourmashSignature) -> Result<SourmashStr> {
    let sig = SourmashSignature::as_rust(ptr);
    Ok(sig.license().into())
}
}

ffi_fn! {
unsafe fn signature_first_mh(ptr: *const SourmashSignature) -> Result<*mut SourmashKmerMinHash> {
    let sig = SourmashSignature::as_rust(ptr);

    if let Some(Sketch::MinHash(mh)) = sig.signatures.get(0) {
        Ok(SourmashKmerMinHash::from_rust(mh.clone()))
    } else {
        // TODO: need to select the correct one
        unimplemented!()
    }
}
}

ffi_fn! {
unsafe fn signature_eq(ptr: *const SourmashSignature, other: *const SourmashSignature) -> Result<bool> {
    let sig = SourmashSignature::as_rust(ptr);
    let other_sig = SourmashSignature::as_rust(other);

    Ok(sig == other_sig)
}
}

ffi_fn! {
unsafe fn signature_save_json(ptr: *const SourmashSignature) -> Result<SourmashStr> {
    let sig = SourmashSignature::as_rust(ptr);
    let st = serde_json::to_string(sig)?;
    Ok(SourmashStr::from_string(st))
}
}

ffi_fn! {
unsafe fn signature_get_mhs(ptr: *const SourmashSignature, size: *mut usize) -> Result<*mut *mut SourmashKmerMinHash> {
    let sig = SourmashSignature::as_rust(ptr);

    let output = sig.sketches();

    // FIXME: how to fit this into the ForeignObject trait?
    let ptr_sigs: Vec<*mut Signature> = output.into_iter().map(|x| {
      Box::into_raw(Box::new(x)) as *mut Signature
    }).collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *mut *mut SourmashKmerMinHash)
}
}

ffi_fn! {
unsafe fn signatures_save_buffer(ptr: *const *const SourmashSignature, size: usize, compression: u8, osize: *mut usize) -> Result<*const u8> {
    // FIXME: review this for ForeignObject

    let sigs = {
        assert!(!ptr.is_null());
        slice::from_raw_parts(ptr, size)
    };

    let rsigs: Vec<&Signature> = sigs.iter().map(|x| SourmashSignature::as_rust(*x)).collect();

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
      serde_json::to_writer(&mut writer, &rsigs)?;
    }

    let b = buffer.into_boxed_slice();
    *osize = b.len();

    Ok(Box::into_raw(b) as *const u8)
}
}

#[no_mangle]
pub unsafe extern "C" fn signature_get_seqlen(ptr: *const SourmashSignature) -> u64 {
    let sig = SourmashSignature::as_rust(ptr);
    sig.seqlen() as u64
}

#[no_mangle]
pub unsafe extern "C" fn signature_set_seqlen(ptr: *mut SourmashSignature, seqlen: u64) {
    let sig = SourmashSignature::as_rust_mut(ptr);
    sig.set_seqlen(seqlen)
}

ffi_fn! {
unsafe fn signatures_load_path(ptr: *const c_char,
                               _ignore_md5sum: bool,
                               ksize: usize,
                               select_moltype: *const c_char,
                               size: *mut usize) -> Result<*mut *mut SourmashSignature> {
    // FIXME use buffer + len instead of cstr
    let buf = {
        assert!(!ptr.is_null());
        CStr::from_ptr(ptr)
    };

    // FIXME take select_moltype as enum
    let moltype: Option<HashFunctions> = if select_moltype.is_null() {
          None
        } else {
          let mol = CStr::from_ptr(select_moltype).to_str()?;
          Some(mol.try_into()?)
    };

    // TODO: implement ignore_md5sum

    let k = match ksize {
      0 => None,
      x => Some(x)
    };

    let (mut input, _) = niffler::from_path(buf.to_str()?)?;
    let filtered_sigs = Signature::load_signatures(&mut input, k, moltype, None)?;

    // FIXME: use the ForeignObject trait, maybe define new method there...
    let ptr_sigs: Vec<*mut SourmashSignature> = filtered_sigs.into_iter().map(|x| {
      Box::into_raw(Box::new(x)) as *mut SourmashSignature
    }).collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *mut *mut SourmashSignature)
}
}

ffi_fn! {
unsafe fn signatures_load_buffer(ptr: *const c_char,
                                 insize: usize,
                                 _ignore_md5sum: bool,
                                 ksize: usize,
                                 select_moltype: *const c_char,
                                 size: *mut usize) -> Result<*mut *mut SourmashSignature> {
    // FIXME use buffer + len instead of cstr
    let buf = {
        assert!(!ptr.is_null());
        slice::from_raw_parts(ptr as *mut u8, insize)
    };

    // FIXME take select_moltype as enum
    let moltype: Option<HashFunctions> = if select_moltype.is_null() {
          None
        } else {
          let mol = CStr::from_ptr(select_moltype).to_str()?;
          Some(mol.try_into()?)
    };

    let k = match ksize {
      0 => None,
      x => Some(x)
    };

    // TODO: implement ignore_md5sum

    let mut reader = io::BufReader::new(buf);
    let filtered_sigs = Signature::load_signatures(&mut reader, k, moltype, None)?;

    // FIXME: use the ForeignObject trait, maybe define new method there...
    let ptr_sigs: Vec<*mut SourmashSignature> = filtered_sigs.into_iter().map(|x| {
      Box::into_raw(Box::new(x)) as *mut SourmashSignature
    }).collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *mut *mut SourmashSignature)
}
}
