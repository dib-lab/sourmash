use std::collections::HashMap;
use std::io::Write;

use crate::index::sbt::{Factory, FromFactory, Node, Update, SBT};
use crate::index::Comparable;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::nodegraph::Nodegraph;
use crate::sketch::Sketch;
use crate::storage::{ReadData, ReadDataError, Storage};
use crate::traits::ToWriter;
use crate::Error;

impl ToWriter for Nodegraph {
    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: Write,
    {
        self.save_to_writer(writer)
    }
}

impl<L: Clone + Default> FromFactory<Node<Nodegraph>> for SBT<Node<Nodegraph>, L> {
    fn factory(&self, name: &str) -> Result<Node<Nodegraph>, Error> {
        match self.factory {
            Factory::GraphFactory { args: (k, t, n) } => {
                let n = Nodegraph::with_tables(t as usize, n as usize, k as usize);

                Ok(Node::builder()
                    .filename(name)
                    .name(name)
                    .metadata(HashMap::default())
                    .storage(self.storage())
                    .data(n)
                    .build())
            }
        }
    }
}

impl Update<Node<Nodegraph>> for Node<Nodegraph> {
    fn update(&self, _other: &mut Node<Nodegraph>) -> Result<(), Error> {
        unimplemented!();
    }
}

impl Update<Node<Nodegraph>> for Signature {
    fn update(&self, parent: &mut Node<Nodegraph>) -> Result<(), Error> {
        // TODO: avoid copy here
        let mut parent_data = parent.data()?.clone();

        if let Sketch::MinHash(sig) = &self.signatures[0] {
            for h in sig.mins() {
                parent_data.count(h);
            }

            let min_n_below = parent
                .metadata
                .entry("min_n_below".into())
                .or_insert(u64::max_value());

            *min_n_below = u64::min(sig.size() as u64, *min_n_below);
            if *min_n_below == 0 {
                *min_n_below = 1
            }
        } else {
            //TODO what if it is not a minhash?
            unimplemented!()
        }

        parent.data = parent_data.into();

        Ok(())
    }
}

impl Comparable<Node<Nodegraph>> for Node<Nodegraph> {
    fn similarity(&self, other: &Node<Nodegraph>) -> f64 {
        let ng: &Nodegraph = self.data().unwrap();
        let ong: &Nodegraph = other.data().unwrap();
        ng.similarity(&ong)
    }

    fn containment(&self, other: &Node<Nodegraph>) -> f64 {
        let ng: &Nodegraph = self.data().unwrap();
        let ong: &Nodegraph = other.data().unwrap();
        ng.containment(&ong)
    }
}

impl Comparable<Signature> for Node<Nodegraph> {
    fn similarity(&self, other: &Signature) -> f64 {
        let ng: &Nodegraph = self.data().unwrap();

        // TODO: select the right signatures...
        if let Sketch::MinHash(sig) = &other.signatures[0] {
            if sig.size() == 0 {
                return 0.0;
            }

            let matches: usize = sig.mins().iter().map(|h| ng.get(*h)).sum();

            let min_n_below = self.metadata["min_n_below"] as f64;

            // This overestimates the similarity, but better than truncating too
            // soon and losing matches
            matches as f64 / min_n_below
        } else {
            //TODO what if it is not a minhash?
            unimplemented!()
        }
    }

    fn containment(&self, other: &Signature) -> f64 {
        let ng: &Nodegraph = self.data().unwrap();

        // TODO: select the right signatures...
        if let Sketch::MinHash(sig) = &other.signatures[0] {
            if sig.size() == 0 {
                return 0.0;
            }

            let matches: usize = sig.mins().iter().map(|h| ng.get(*h)).sum();

            matches as f64 / sig.size() as f64
        } else {
            //TODO what if it is not a minhash?
            unimplemented!()
        }
    }
}

impl ReadData<Nodegraph> for Node<Nodegraph> {
    fn data(&self) -> Result<&Nodegraph, Error> {
        if let Some(storage) = &self.storage {
            Ok(self.data.get_or_init(|| {
                let raw = storage.load(&self.filename).unwrap();
                Nodegraph::from_reader(&mut &raw[..]).unwrap()
            }))
        } else if let Some(data) = self.data.get() {
            Ok(data)
        } else {
            Err(ReadDataError::LoadError.into())
        }
    }
}

#[cfg(test)]
mod test {
    use std::convert::TryInto;
    use std::fs::File;
    use std::io::{BufReader, Seek, SeekFrom};
    use std::path::PathBuf;

    use assert_matches::assert_matches;

    use super::Factory;

    use crate::index::linear::LinearIndex;
    use crate::index::sbt::scaffold;
    use crate::index::search::{search_minhashes, search_minhashes_containment};
    use crate::index::{Index, SigStore, MHBT};
    use crate::signature::Signature;
    use crate::storage::ReadData;

    #[test]
    fn save_sbt() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/v5.sbt.json");

        let mut sbt = MHBT::from_path(filename).expect("Loading error");

        let mut tmpfile = tempfile::NamedTempFile::new().unwrap();
        sbt.save_file(tmpfile.path(), None).unwrap();

        tmpfile.seek(SeekFrom::Start(0)).unwrap();
    }

    #[test]
    fn load_sbt() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/v5.sbt.json");

        let sbt = MHBT::from_path(filename).expect("Loading error");

        assert_eq!(sbt.d, 2);
        //assert_eq!(sbt.storage.backend, "FSStorage");
        //assert_eq!(sbt.storage.args["path"], ".sbt.v5");
        //assert_matches!(&sbt.storage, <dyn Storage as Trait>::FSStorage(args) => {
        //    assert_eq!(args, &[1, 100000, 4]);
        //});
        assert_matches!(&sbt.factory, Factory::GraphFactory { args } => {
            assert_eq!(args, &(1, 100000.0, 4));
        });

        println!("sbt leaves {:?} {:?}", sbt.leaves.len(), sbt.leaves);

        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/.sbt.v3/60f7e23c24a8d94791cc7a8680c493f9");

        let mut reader = BufReader::new(File::open(filename).unwrap());
        let sigs = Signature::load_signatures(
            &mut reader,
            Some(31),
            Some("DNA".try_into().unwrap()),
            None,
        )
        .unwrap();
        let leaf = sigs[0].clone();

        let results = sbt.find(search_minhashes, &leaf, 0.5).unwrap();
        assert_eq!(results.len(), 1);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let results = sbt.find(search_minhashes, &leaf, 0.1).unwrap();
        assert_eq!(results.len(), 2);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let mut linear = LinearIndex::builder().storage(sbt.storage()).build();
        for l in &sbt.leaves {
            linear.insert(l.1.data().unwrap().clone()).unwrap();
        }

        let datasets = linear.signatures();
        println!("linear leaves {:?} {:?}", datasets.len(), datasets);

        let results = linear.find(search_minhashes, &leaf, 0.5).unwrap();
        assert_eq!(results.len(), 1);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let results = linear.find(search_minhashes, &leaf, 0.1).unwrap();
        assert_eq!(results.len(), 2);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let results = linear
            .find(search_minhashes_containment, &leaf, 0.5)
            .unwrap();
        assert_eq!(results.len(), 2);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let results = linear
            .find(search_minhashes_containment, &leaf, 0.1)
            .unwrap();
        assert_eq!(results.len(), 4);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);
    }

    #[test]
    #[ignore]
    fn roundtrip_sbt() -> Result<(), Box<dyn std::error::Error>> {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/v5.sbt.json");

        let sbt = MHBT::from_path(filename)?;

        assert_eq!(sbt.d, 2);
        //assert_eq!(sbt.storage.backend, "FSStorage");
        //assert_eq!(sbt.storage.args["path"], ".sbt.v5");
        //assert_matches!(&sbt.storage, <dyn Storage as Trait>::FSStorage(args) => {
        //    assert_eq!(args, &[1, 100000, 4]);
        //});
        assert_matches!(&sbt.factory, Factory::GraphFactory { args } => {
            assert_eq!(args, &(1, 100000.0, 4));
        });

        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/.sbt.v3/60f7e23c24a8d94791cc7a8680c493f9");

        let mut reader = BufReader::new(File::open(filename)?);
        let sigs = Signature::load_signatures(
            &mut reader,
            Some(31),
            Some("DNA".try_into().unwrap()),
            None,
        )?;
        let sig_data = sigs[0].clone();

        let leaf: SigStore<_> = sig_data.into();

        let results = sbt.find(search_minhashes, &leaf, 0.5)?;
        assert_eq!(results.len(), 1);
        //println!("results: {:?}", results);
        //println!("leaf: {:?}", leaf);

        let results = sbt.find(search_minhashes, &leaf, 0.1)?;
        assert_eq!(results.len(), 2);
        //println!("results: {:?}", results);
        //println!("leaf: {:?}", leaf);

        println!("sbt internal {:?} {:?}", sbt.nodes.len(), sbt.nodes);
        println!("sbt leaves {:?} {:?}", sbt.leaves.len(), sbt.leaves);

        let mut new_sbt: MHBT = MHBT::builder().storage(None).build();
        let datasets = sbt.signatures();
        for l in datasets {
            new_sbt.insert(l)?;
        }

        for (i, node) in &sbt.nodes {
            assert_eq!(node.data().unwrap(), new_sbt.nodes[i].data().unwrap());
        }

        assert_eq!(new_sbt.signature_refs().len(), 7);
        println!("new_sbt internal {:?} {:?}", sbt.nodes.len(), sbt.nodes);
        println!("new_sbt leaves {:?} {:?}", sbt.leaves.len(), sbt.leaves);

        let results = new_sbt.find(search_minhashes, &leaf, 0.5)?;
        //println!("results: {:?}", results);
        //println!("leaf: {:?}", leaf);
        assert_eq!(results.len(), 1);

        let results = new_sbt.find(search_minhashes, &leaf, 0.1)?;
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);
        assert_eq!(results.len(), 2);

        let results = new_sbt.find(search_minhashes_containment, &leaf, 0.5)?;
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);
        assert_eq!(results.len(), 2);

        let results = new_sbt.find(search_minhashes_containment, &leaf, 0.1)?;
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);
        assert_eq!(results.len(), 4);

        Ok(())
    }

    #[test]
    fn scaffold_sbt() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/v5.sbt.json");

        let sbt = MHBT::from_path(filename).expect("Loading error");

        let new_sbt: MHBT = scaffold(sbt.leaves(), sbt.storage());

        assert_eq!(new_sbt.signatures().len(), 7);
    }

    #[test]
    fn load_v4() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/v4.sbt.json");

        let _sbt = MHBT::from_path(filename).expect("Loading error");
    }

    #[test]
    fn load_v5() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/v5.sbt.json");

        let _sbt = MHBT::from_path(filename).expect("Loading error");
    }
}
