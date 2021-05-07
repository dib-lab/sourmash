use crate::ffi::utils::ForeignObject;
use crate::index::storage::MemStorage;

pub struct SourmashMemStorage;

impl ForeignObject for SourmashMemStorage {
    type RustObject = MemStorage;
}

#[no_mangle]
pub unsafe extern "C" fn memstorage_new() -> *mut SourmashMemStorage {
    let memstorage = MemStorage::default();

    SourmashMemStorage::from_rust(memstorage)
}

#[no_mangle]
pub unsafe extern "C" fn memstorage_free(ptr: *mut SourmashMemStorage) {
    SourmashMemStorage::drop(ptr);
}
