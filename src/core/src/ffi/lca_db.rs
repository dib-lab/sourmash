use std::ffi::CStr;
use std::os::raw::c_char;
use std::path::Path;

use crate::ffi::utils::{ForeignObject, SourmashStr};
use crate::index::lca_db::{LcaDB, Lineage, lineage_to_vec, zip_lineage};
use crate::signature::Signature;
use crate::ffi::signature::SourmashSignature;
use std::collections::{HashMap, BTreeMap};


pub struct SourmashLcaDatabase;

impl ForeignObject for SourmashLcaDatabase {
    type RustObject = LcaDB;
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_new_with_params(
    ksize: u32,
    scaled: u64,
    filename_char: *const c_char,
    moltype_char: *const c_char,
) -> *mut SourmashLcaDatabase {
    let filename = {
        assert!(!filename_char.is_null());
        CStr::from_ptr(filename_char).to_str().unwrap()
    };
    let moltype = {
        assert!(!moltype_char.is_null());
        CStr::from_ptr(moltype_char).to_str().unwrap()
    };
    
    let lca_db = LcaDB::new_with_params(Some(ksize), Some(scaled), Some(filename), Some(moltype));

    SourmashLcaDatabase::from_rust(lca_db)
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_ksize(ptr: *const SourmashLcaDatabase) -> u32 {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    lca_db.ksize()
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_scaled(ptr: *const SourmashLcaDatabase) -> u64 {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    lca_db.scaled()
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_moltype(ptr: *const SourmashLcaDatabase) -> SourmashStr {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    lca_db.moltype().into()
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_filename(ptr: *const SourmashLcaDatabase) -> SourmashStr {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    lca_db.filename().into()
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_ident_to_name(ptr: *const SourmashLcaDatabase) -> SourmashStr {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    let result = serde_json::to_string(&lca_db.ident_to_name()).unwrap();

    result.into()
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_ident_to_idx(ptr: *const SourmashLcaDatabase) -> SourmashStr {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    let result = serde_json::to_string(&lca_db.ident_to_idx()).unwrap();

    result.into()
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_idx_to_lid(ptr: *const SourmashLcaDatabase) -> SourmashStr {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    let result = serde_json::to_string(&lca_db.idx_to_lid()).unwrap();

    result.into()
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_lineage_to_lid(ptr: *const SourmashLcaDatabase) -> SourmashStr {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);

    let mut lineage_to_lid_vec: HashMap<Vec<Vec<String>>, u32> = HashMap::new();
    for (lineage, lid) in lca_db.lineage_to_lid() {
        lineage_to_lid_vec.insert(lineage_to_vec(lineage), lid);
    }

    let result = serde_json::to_string(&lineage_to_lid_vec).unwrap();

    result.into()
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_lid_to_lineage(ptr: *const SourmashLcaDatabase) -> SourmashStr {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);

    let mut lid_to_lineage_vec: HashMap<u32, Vec<Vec<String>>> = HashMap::new();
    for (lid, lineage) in lca_db.lid_to_lineage() {
        lid_to_lineage_vec.insert(lid, lineage_to_vec(lineage));
    }

    let result = serde_json::to_string(&lid_to_lineage_vec).unwrap();

    result.into()
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_hashval_to_idx(ptr: *const SourmashLcaDatabase) -> SourmashStr {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);

    let result = serde_json::to_string(&lca_db.hashval_to_idx()).unwrap();

    result.into()
}

ffi_fn! {
unsafe fn lcadb_get_lineage_from_idx(
    ptr: *mut SourmashLcaDatabase, 
    idx: u32, 
) -> Result<SourmashStr> {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);

    let lineage: Lineage = match lca_db.idx_to_lid().get(&idx) {
        Some(lid) => lca_db.lid_to_lineage()[&lid].clone(),
        None => BTreeMap::new(),
    };
  
    // let result = serde_json::to_string(&lineage_to_vec(lineage)).unwrap();
    let result = zip_lineage(&lineage).unwrap_or(serde_json::to_string(&lineage).unwrap());

    Ok(result.into())
}
}

ffi_fn! {
unsafe fn lcadb_get_match_size(
    ptr: *mut SourmashLcaDatabase, 
    best_idx: u32,
) -> Result<u32> {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);

    let mut match_size = 0;
    for (_hashval, idx) in lca_db.hashval_to_idx() {
        if idx.contains(&best_idx) {
            match_size += 1;
        }
    }

    Ok(match_size)
}
}

ffi_fn! {
unsafe fn lcadb_best_name(
    ptr: *mut SourmashLcaDatabase, 
    best_idx: u32,
) -> Result<SourmashStr> {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);

    let mut name = String::new();
    if let Some(ident) = lca_db.idx_to_ident().get(&best_idx) {
        name = lca_db.ident_to_name()[&ident.clone()].clone();
    }

    Ok(name.into())
}
}

ffi_fn! {
unsafe fn lcadb_get_lineage_assignments(
    ptr: *mut SourmashLcaDatabase, 
    hashval: u64, 
    size: *mut usize,
) -> Result<*const SourmashStr> {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);

    let assignments = lca_db.get_lineage_assignments(hashval).unwrap();

    let x: Vec<SourmashStr> = assignments.into_iter().map(|lineage| {
        zip_lineage(&lineage).unwrap_or("".to_string()).into()
    }).collect();

    let b = x.into_boxed_slice();
    *size = b.len();
    let result = Box::into_raw(b) as *const SourmashStr;

    Ok(result)
}
}

ffi_fn! {
unsafe fn lcadb_downsample_scaled(
    ptr: *mut SourmashLcaDatabase, 
    new_scaled: u64,
) -> Result<()> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    
    lca_db.downsample_scaled(new_scaled)?;

    Ok(())
}
}

ffi_fn! {
unsafe fn lcadb_insert(
    ptr: *mut SourmashLcaDatabase,
    sig_ptr:  *const SourmashSignature,
    ident_char: *const c_char,
    lineage_char: *const c_char,
) -> Result<u32> {

    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    let sig = SourmashSignature::as_rust(sig_ptr);
    let ident = {
        assert!(!ident_char.is_null());
        CStr::from_ptr(ident_char)
    };

    //set lineage
    let data = CStr::from_ptr(lineage_char).to_str().unwrap();
    let lineage = if data != "" {
        let lineage_temp: Lineage = serde_json::from_str(data)?;
        Some(lineage_temp)
    } else {
        None
    };
    
    lca_db.insert(sig, Some(ident.to_str().unwrap()), lineage)
}
}

// WORKS!
ffi_fn! {
unsafe fn lcadb_signatures(
    ptr: *mut SourmashLcaDatabase,
    size: *mut usize,
) -> Result<*mut *mut SourmashSignature> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);

    // signatures() -> Result<Vec<Signature>>
    let sigs = lca_db.signatures()?;

    // FIXME: use the ForeignObject trait, maybe define new method there...
    let ptr_sigs: Vec<*mut SourmashSignature> = sigs.into_iter().map(|x| {
      Box::into_raw(Box::new(x)) as *mut SourmashSignature
    }).collect();

    
    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *mut *mut SourmashSignature)
}
}

ffi_fn! {
unsafe fn lcadb_load(
    filepath_char: *const c_char,
) -> Result<*mut SourmashLcaDatabase> {

    let buf = {
        assert!(!filepath_char.is_null());
        CStr::from_ptr(filepath_char)
    };
    let filename = &buf.to_str()?;
    
    let path = Path::new(filename);
    let lca_db = LcaDB::load(path, filename).unwrap();

    Ok(SourmashLcaDatabase::from_rust(lca_db))
}
}

ffi_fn! {
unsafe fn lcadb_load_db(
    filepath_char: *const c_char,
) -> Result<*mut SourmashLcaDatabase> {

    let buf = {
        assert!(!filepath_char.is_null());
        CStr::from_ptr(filepath_char)
    };
    let filename = &buf.to_str()?;
    
    let (mut input, _) = niffler::from_path(filename)?;
    dbg!(&filename);
    let lca_db = LcaDB::load_db(input, filename).unwrap();

    Ok(SourmashLcaDatabase::from_rust(lca_db))
}
}

ffi_fn! {
unsafe fn lcadb_save(
    ptr: *mut SourmashLcaDatabase, 
    filepath_char: *const c_char,
) -> Result<()> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    let filepath = {
        assert!(!filepath_char.is_null());
        CStr::from_ptr(filepath_char)
    };
    let path = Path::new(filepath.to_str().unwrap());
    lca_db.save(path)
}
}

#[repr(C)]
pub struct FFISearchResults {
    score: f32,
    sig: *mut SourmashSignature,
    filename: SourmashStr,
}

impl FFISearchResults {
    pub fn new(score: f32, sig: Signature, filename: String) -> FFISearchResults {
        let filename: SourmashStr  = filename.into();
        let sig = unsafe { SourmashSignature::from_rust(sig) };

        FFISearchResults {
            score,
            sig,
            filename,
        }
    }
}

ffi_fn! {
unsafe fn lcadb_search(
    ptr: *mut SourmashLcaDatabase,
    query: *const SourmashSignature,
    threshold: f32,
    do_containment: bool,
    ignore_abundance: bool,
    osize: *mut usize
) -> Result<*mut FFISearchResults> {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    let sig = SourmashSignature::as_rust(query);

    let searchresults: Vec<(f32, Signature, String)> = lca_db.search(sig.clone(), threshold, do_containment, ignore_abundance)?;
    
    let result_vec: Vec<FFISearchResults> = searchresults.into_iter().map(|x| FFISearchResults::new(x.0, x.1, x.2)).collect();
    
    let b = result_vec.into_boxed_slice();
    *osize = b.len();

    Ok(Box::into_raw(b) as *mut FFISearchResults)
}
}

ffi_fn! {
unsafe fn lcadb_gather(
    ptr: *mut SourmashLcaDatabase,
    query: *const SourmashSignature,
    threshold_bp: f32,
    osize: *mut usize
) -> Result<*mut FFISearchResults> {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    let sig = SourmashSignature::as_rust(query);

    let searchresults: Vec<(f32, Signature, String)> = lca_db.gather(sig.clone(), threshold_bp).unwrap();

    let result_vec: Vec<FFISearchResults> = searchresults.into_iter().map(|x| FFISearchResults::new(x.0, x.1, x.2)).collect();
    
    let b = result_vec.into_boxed_slice();
    *osize = b.len();

    Ok(Box::into_raw(b) as *mut FFISearchResults)
}
}
