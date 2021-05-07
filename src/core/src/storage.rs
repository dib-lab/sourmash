use std::fs::{DirBuilder, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use serde::{Deserialize, Serialize};
use thiserror::Error;
use typed_builder::TypedBuilder;

use crate::Error;

#[derive(Clone)]
pub struct InnerStorage(Arc<Mutex<dyn Storage>>);

impl InnerStorage {
    pub fn new(inner: impl Storage + 'static) -> InnerStorage {
        InnerStorage(Arc::new(Mutex::new(inner)))
    }
}

impl Storage for InnerStorage {
    fn save(&self, path: &str, content: &[u8]) -> Result<String, Error> {
        self.0.save(path, content)
    }
    fn load(&self, path: &str) -> Result<Vec<u8>, Error> {
        self.0.load(path)
    }
    fn args(&self) -> StorageArgs {
        self.0.args()
    }
}

#[derive(Debug, Error)]
pub enum StorageError {
    #[error("Path can't be empty")]
    EmptyPathError,
}

#[derive(Debug, Error)]
pub enum ReadDataError {
    #[error("Could not load data")]
    LoadError,
}

/// Implemented by anything that wants to read specific data from a storage.
pub trait ReadData<D> {
    fn data(&self) -> Result<&D, Error>;
}

#[derive(Serialize, Deserialize)]
pub(crate) struct StorageInfo {
    pub backend: String,
    pub args: StorageArgs,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum StorageArgs {
    FSStorage { path: String },
}

impl From<&StorageArgs> for FSStorage {
    fn from(other: &StorageArgs) -> FSStorage {
        match other {
            StorageArgs::FSStorage { path } => {
                let mut fullpath = PathBuf::new();
                fullpath.push(".");
                fullpath.push(path);

                FSStorage {
                    fullpath,
                    subdir: path.clone(),
                }
            }
        }
    }
}

/// An abstraction for any place where we can store data.
pub trait Storage {
    /// Save bytes into path
    fn save(&self, path: &str, content: &[u8]) -> Result<String, Error>;

    /// Load bytes from path
    fn load(&self, path: &str) -> Result<Vec<u8>, Error>;

    /// Args for initializing a new Storage
    fn args(&self) -> StorageArgs;
}

impl<L> Storage for Mutex<L>
where
    L: ?Sized + Storage,
{
    fn save(&self, path: &str, content: &[u8]) -> Result<String, Error> {
        self.lock().unwrap().save(path, content)
    }

    fn load(&self, path: &str) -> Result<Vec<u8>, Error> {
        self.lock().unwrap().load(path)
    }

    fn args(&self) -> StorageArgs {
        self.lock().unwrap().args()
    }
}

/// Store files locally into a directory
#[derive(TypedBuilder, Debug, Clone, Default)]
pub struct FSStorage {
    /// absolute path for the directory where data is saved.
    fullpath: PathBuf,
    subdir: String,
}

impl FSStorage {
    pub fn new(location: &str, subdir: &str) -> FSStorage {
        let mut fullpath = PathBuf::new();
        fullpath.push(location);
        fullpath.push(subdir);

        FSStorage {
            fullpath,
            subdir: subdir.into(),
        }
    }

    pub fn set_base(&mut self, location: &str) {
        let mut fullpath = PathBuf::new();
        fullpath.push(location);
        fullpath.push(&self.subdir);
        self.fullpath = fullpath;
    }
}

impl Storage for FSStorage {
    fn save(&self, path: &str, content: &[u8]) -> Result<String, Error> {
        if path.is_empty() {
            return Err(StorageError::EmptyPathError.into());
        }

        let fpath = self.fullpath.join(path);
        DirBuilder::new()
            .recursive(true)
            .create(fpath.parent().unwrap())?;

        let file = File::create(&fpath)?;
        let mut buf_writer = BufWriter::new(file);
        buf_writer.write_all(content)?;
        Ok(path.into())
    }

    fn load(&self, path: &str) -> Result<Vec<u8>, Error> {
        let path = self.fullpath.join(path);
        let file = File::open(path)?;
        let mut buf_reader = BufReader::new(file);
        let mut contents = Vec::new();
        buf_reader.read_to_end(&mut contents)?;
        Ok(contents)
    }

    fn args(&self) -> StorageArgs {
        StorageArgs::FSStorage {
            path: self.subdir.clone(),
        }
    }
}
