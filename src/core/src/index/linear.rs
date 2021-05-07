use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;
use std::path::PathBuf;

use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

use crate::index::{Comparable, DatasetInfo, Index, SigStore};
use crate::storage::{FSStorage, InnerStorage, ReadData, Storage, StorageInfo};
use crate::traits::ToWriter;
use crate::Error;

#[derive(TypedBuilder)]
pub struct LinearIndex<L> {
    #[builder(default)]
    storage: Option<InnerStorage>,

    #[builder(default)]
    datasets: Vec<SigStore<L>>,
}

#[derive(Serialize, Deserialize)]
struct LinearInfo<L> {
    version: u32,
    storage: StorageInfo,
    leaves: Vec<L>,
}

impl<'a, L> Index<'a> for LinearIndex<L>
where
    L: Clone + Comparable<L> + 'a,
    SigStore<L>: From<L>,
{
    type Item = L;
    //type SignatureIterator = std::slice::Iter<'a, Self::Item>;

    fn insert(&mut self, node: L) -> Result<(), Error> {
        self.datasets.push(node.into());
        Ok(())
    }

    fn save<P: AsRef<Path>>(&self, _path: P) -> Result<(), Error> {
        /*
        let file = File::create(path)?;
        match serde_json::to_writer(file, &self) {
            Ok(_) => Ok(()),
            Err(_) => Err(SourmashError::SerdeError.into()),
        }
        */
        unimplemented!()
    }

    fn load<P: AsRef<Path>>(_path: P) -> Result<(), Error> {
        unimplemented!()
    }

    fn signatures(&self) -> Vec<Self::Item> {
        self.datasets
            .iter()
            .map(|x| x.data.get().unwrap().clone())
            .collect()
    }

    fn signature_refs(&self) -> Vec<&Self::Item> {
        self.datasets
            .iter()
            .map(|x| x.data.get().unwrap())
            .collect()
    }

    /*
    fn iter_signatures(&'a self) -> Self::SignatureIterator {
        self.datasets.iter()
    }
    */
}

impl<L> LinearIndex<L>
where
    L: ToWriter,
    SigStore<L>: ReadData<L>,
{
    pub fn save_file<P: AsRef<Path>>(
        &mut self,
        path: P,
        storage: Option<InnerStorage>,
    ) -> Result<(), Error> {
        let ref_path = path.as_ref();
        let mut basename = ref_path.file_name().unwrap().to_str().unwrap().to_owned();
        if basename.ends_with(".sbt.json") {
            basename = basename.replace(".sbt.json", "");
        }
        let location = ref_path.parent().unwrap();

        let storage = match storage {
            Some(s) => s,
            None => {
                let subdir = format!(".linear.{}", basename);
                InnerStorage::new(FSStorage::new(location.to_str().unwrap(), &subdir))
            }
        };

        let args = storage.args();
        let storage_info = StorageInfo {
            backend: "FSStorage".into(),
            args,
        };

        let info: LinearInfo<DatasetInfo> = LinearInfo {
            storage: storage_info,
            version: 5,
            leaves: self
                .datasets
                .iter_mut()
                .map(|l| {
                    // Trigger data loading
                    let _: &L = (*l).data().unwrap();

                    // set storage to new one
                    l.storage = Some(storage.clone());

                    let filename = (*l).save(&l.filename).unwrap();

                    DatasetInfo {
                        filename,
                        name: l.name.clone(),
                        metadata: l.metadata.clone(),
                    }
                })
                .collect(),
        };

        let file = File::create(path)?;
        serde_json::to_writer(file, &info)?;

        Ok(())
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<LinearIndex<L>, Error> {
        let file = File::open(&path)?;
        let mut reader = BufReader::new(file);

        // TODO: match with available Storage while we don't
        // add a function to build a Storage from a StorageInfo
        let mut basepath = PathBuf::new();
        basepath.push(path);
        basepath.canonicalize()?;

        let linear = LinearIndex::<L>::from_reader(&mut reader, &basepath.parent().unwrap())?;
        Ok(linear)
    }

    pub fn from_reader<R, P>(rdr: R, path: P) -> Result<LinearIndex<L>, Error>
    where
        R: Read,
        P: AsRef<Path>,
    {
        // TODO: check https://serde.rs/enum-representations.html for a
        // solution for loading v4 and v5
        let linear: LinearInfo<DatasetInfo> = serde_json::from_reader(rdr)?;

        // TODO: support other storages
        let mut st: FSStorage = (&linear.storage.args).into();
        st.set_base(path.as_ref().to_str().unwrap());
        let storage = InnerStorage::new(st);

        Ok(LinearIndex {
            storage: Some(storage.clone()),
            datasets: linear
                .leaves
                .into_iter()
                .map(|l| {
                    let mut v: SigStore<L> = l.into();
                    v.storage = Some(storage.clone());
                    v
                })
                .collect(),
        })
    }

    pub fn storage(&self) -> Option<InnerStorage> {
        self.storage.clone()
    }
}
