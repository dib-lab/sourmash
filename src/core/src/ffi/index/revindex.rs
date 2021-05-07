use std::path::PathBuf;
use std::slice;

use crate::index::greyhound::RevIndex;
use crate::index::Index;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::minhash::KmerMinHash;
use crate::sketch::Sketch;

use crate::ffi::index::SourmashSearchResult;
use crate::ffi::minhash::SourmashKmerMinHash;
use crate::ffi::signature::SourmashSignature;
use crate::ffi::utils::{ForeignObject, SourmashStr};

pub struct SourmashRevIndex;

impl ForeignObject for SourmashRevIndex {
    type RustObject = RevIndex;
}

ffi_fn! {
  unsafe fn revindex_new_with_paths(
      search_sigs_ptr: *const *const SourmashStr,
      insigs: usize,
      template_ptr: *const SourmashKmerMinHash,
      threshold: usize,
      queries_ptr: *const *const SourmashKmerMinHash,
      inqueries: usize,
      keep_sigs: bool,
  ) -> Result<*mut SourmashRevIndex> {
    let search_sigs: Vec<PathBuf> = {
      assert!(!search_sigs_ptr.is_null());
        slice::from_raw_parts(search_sigs_ptr, insigs).iter().map(|path| {
          let mut new_path = PathBuf::new();
          new_path.push(SourmashStr::as_rust(*path).as_str());
          new_path}
          ).collect()
    };

    let template = {
      assert!(!template_ptr.is_null());
      //TODO: avoid clone here
      Sketch::MinHash(SourmashKmerMinHash::as_rust(template_ptr).clone())
    };

    let queries_vec: Vec<KmerMinHash>;
    let queries: Option<&[KmerMinHash]> = if queries_ptr.is_null() {
      None
    } else {
        queries_vec =
          slice::from_raw_parts(queries_ptr, inqueries).iter().map(|mh_ptr|
            // TODO: avoid this clone
          SourmashKmerMinHash::as_rust(*mh_ptr).clone()).collect();
        Some(queries_vec.as_ref())
    };
      let revindex = RevIndex::new(
        search_sigs.as_ref(),
        &template,
        threshold,
        queries,
        keep_sigs
      );
      Ok(SourmashRevIndex::from_rust(revindex))
  }
}

ffi_fn! {
  unsafe fn revindex_new_with_sigs(
      search_sigs_ptr: *const *const SourmashSignature,
      insigs: usize,
      template_ptr: *const SourmashKmerMinHash,
      threshold: usize,
      queries_ptr: *const *const SourmashKmerMinHash,
      inqueries: usize,
  ) -> Result<*mut SourmashRevIndex> {
    let search_sigs: Vec<Signature> = {
      assert!(!search_sigs_ptr.is_null());
        slice::from_raw_parts(search_sigs_ptr, insigs).iter().map(|sig|
          SourmashSignature::as_rust(*sig)
          ).cloned().collect()
    };

    let template = {
      assert!(!template_ptr.is_null());
      //TODO: avoid clone here
      Sketch::MinHash(SourmashKmerMinHash::as_rust(template_ptr).clone())
    };

    let queries_vec: Vec<KmerMinHash>;
    let queries: Option<&[KmerMinHash]> = if queries_ptr.is_null() {
      None
    } else {
        queries_vec =
          slice::from_raw_parts(queries_ptr, inqueries).iter().map(|mh_ptr|
            // TODO: avoid this clone
          SourmashKmerMinHash::as_rust(*mh_ptr).clone()).collect();
        Some(queries_vec.as_ref())
    };
      let revindex = RevIndex::new_with_sigs(
        search_sigs,
        &template,
        threshold,
        queries,
      );
      Ok(SourmashRevIndex::from_rust(revindex))
  }
}

#[no_mangle]
pub unsafe extern "C" fn revindex_free(ptr: *mut SourmashRevIndex) {
    SourmashRevIndex::drop(ptr);
}

ffi_fn! {
unsafe fn revindex_search(
    ptr: *const SourmashRevIndex,
    sig_ptr: *const SourmashSignature,
    threshold: f64,
    do_containment: bool,
    _ignore_abundance: bool,
    size: *mut usize
) -> Result<*const *const SourmashSearchResult> {
    let revindex = SourmashRevIndex::as_rust(ptr);
    let sig = SourmashSignature::as_rust(sig_ptr);

    if sig.signatures.is_empty() {
        *size = 0;
        return Ok(std::ptr::null::<*const SourmashSearchResult>());
    }

    let mh = if let Sketch::MinHash(mh) = &sig.signatures[0] {
        mh
    } else {
        // TODO: what if it is not a mh?
        unimplemented!()
    };

    let results: Vec<(f64, Signature, String)> = revindex
        .find_signatures(mh, threshold, do_containment, true)?
        .into_iter()
        .collect();

    // FIXME: use the ForeignObject trait, maybe define new method there...
    let ptr_sigs: Vec<*const SourmashSearchResult> = results.into_iter().map(|x| {
      Box::into_raw(Box::new(x)) as *const SourmashSearchResult
    }).collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *const *const SourmashSearchResult)
}
}

ffi_fn! {
unsafe fn revindex_gather(
    ptr: *const SourmashRevIndex,
    sig_ptr: *const SourmashSignature,
    threshold: f64,
    _do_containment: bool,
    _ignore_abundance: bool,
    size: *mut usize
) -> Result<*const *const SourmashSearchResult> {
    let revindex = SourmashRevIndex::as_rust(ptr);
    let sig = SourmashSignature::as_rust(sig_ptr);

    if sig.signatures.is_empty() {
        *size = 0;
        return Ok(std::ptr::null::<*const SourmashSearchResult>());
    }

    let mh = if let Sketch::MinHash(mh) = &sig.signatures[0] {
        mh
    } else {
        // TODO: what if it is not a mh?
        unimplemented!()
    };

    // TODO: proper threshold calculation
    let threshold: usize = (threshold * (mh.size() as f64)) as _;

    let counter = revindex.counter_for_query(&mh);
        dbg!(&counter);

    let results: Vec<(f64, Signature, String)> = revindex
        .gather(counter, threshold, mh)
        .unwrap() // TODO: proper error handling
        .into_iter()
        .map(|r| {
            let filename = r.filename().to_owned();
            let sig = r.get_match();
            (r.f_match(), sig, filename)
        })
        .collect();


    // FIXME: use the ForeignObject trait, maybe define new method there...
    let ptr_sigs: Vec<*const SourmashSearchResult> = results.into_iter().map(|x| {
      Box::into_raw(Box::new(x)) as *const SourmashSearchResult
    }).collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *const *const SourmashSearchResult)
}
}

#[no_mangle]
pub unsafe extern "C" fn revindex_scaled(ptr: *const SourmashRevIndex) -> u64 {
    let revindex = SourmashRevIndex::as_rust(ptr);
    if let Sketch::MinHash(mh) = revindex.template() {
        mh.scaled()
    } else {
        unimplemented!()
    }
}

ffi_fn! {
unsafe fn revindex_signatures(ptr: *const SourmashRevIndex,
                              size: *mut usize) -> Result<*mut *mut SourmashSignature> {
    let revindex = SourmashRevIndex::as_rust(ptr);

    let sigs = revindex.signatures();

    // FIXME: use the ForeignObject trait, maybe define new method there...
    let ptr_sigs: Vec<*mut SourmashSignature> = sigs.into_iter().map(|x| {
      Box::into_raw(Box::new(x)) as *mut SourmashSignature
    }).collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *mut *mut SourmashSignature)
}
}
