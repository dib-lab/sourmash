use std::collections::HashMap;
use std::io::{Read, Write};

use failure::Error;
use mqf::MQF;
use tempfile;

use crate::index::sbt::{Factory, FromFactory, Node, Update, SBT};
use crate::index::storage::{ReadData, ReadDataError, ToWriter};
use crate::index::Comparable;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::Sketch;

impl ToWriter for MQF {
    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: Write,
    {
        // TODO: using tempfile for now, but ideally want to avoid that
        let mut tmpfile = tempfile::NamedTempFile::new()?;
        self.serialize(tmpfile.path()).unwrap(); // TODO: convert this to a proper error

        let mut buffer = Vec::new();
        tmpfile.read_to_end(&mut buffer)?;
        writer.write_all(&buffer)?;

        Ok(())
    }
}

impl ReadData<MQF> for Node<MQF> {
    fn data(&self) -> Result<&MQF, Error> {
        if let Some(storage) = &self.storage {
            Ok(self.data.get_or_init(|| {
                let raw = storage.load(&self.filename).unwrap();

                // TODO: using tempfile for now, but ideally want to avoid that
                let mut tmpfile = tempfile::NamedTempFile::new().unwrap();
                tmpfile.write_all(&raw[..]).unwrap();

                MQF::deserialize(tmpfile.path()).unwrap()
            }))
        } else if let Some(data) = self.data.get() {
            Ok(data)
        } else {
            Err(ReadDataError::LoadError.into())
        }
    }
}

impl<L: Clone + Default> FromFactory<Node<MQF>> for SBT<Node<MQF>, L> {
    fn factory(&self, name: &str) -> Result<Node<MQF>, Error> {
        match self.factory {
            Factory::GraphFactory { .. } => {
                // TODO: define a new factory and use it properly
                let n = MQF::new(1, 18);

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

impl Update<Node<MQF>> for Node<MQF> {
    fn update(&self, _other: &mut Node<MQF>) -> Result<(), Error> {
        unimplemented!();
    }
}

impl Update<Node<MQF>> for Signature {
    fn update(&self, parent: &mut Node<MQF>) -> Result<(), Error> {
        // TODO: avoid copy here
        let mut parent_data = parent.data()?.clone();

        if let Sketch::MinHash(sig) = &self.signatures[0] {
            sig.mins
                .iter()
                .for_each(|h| parent_data.insert(*h % u64::pow(2, 26), 1));

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

impl Comparable<Node<MQF>> for Node<MQF> {
    fn similarity(&self, other: &Node<MQF>) -> f64 {
        let _ng: &MQF = self.data().unwrap();
        let _ong: &MQF = other.data().unwrap();
        unimplemented!();
        //ng.similarity(&ong)
    }

    fn containment(&self, other: &Node<MQF>) -> f64 {
        let _ng: &MQF = self.data().unwrap();
        let _ong: &MQF = other.data().unwrap();
        unimplemented!();
        //ng.containment(&ong)
    }
}

impl Comparable<Signature> for Node<MQF> {
    fn similarity(&self, other: &Signature) -> f64 {
        let ng: &MQF = self.data().unwrap();

        // TODO: select the right signatures...
        if let Sketch::MinHash(sig) = &other.signatures[0] {
            if sig.size() == 0 {
                return 0.0;
            }

            let matches: usize = sig
                .mins
                .iter()
                .filter(|h| ng.count_key(**h % u64::pow(2, 26)) > 0)
                //.filter(|h| ng.count_key(**h) > 0)
                .count();

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
        let ng: &MQF = self.data().unwrap();

        // TODO: select the right signatures...
        if let Sketch::MinHash(sig) = &other.signatures[0] {
            if sig.size() == 0 {
                return 0.0;
            }

            let matches: usize = sig
                .mins
                .iter()
                .filter(|h| ng.count_key(**h % u64::pow(2, 26)) > 0)
                //.filter(|h| ng.count_key(**h) > 0)
                .count();

            matches as f64 / sig.size() as f64
        } else {
            //TODO what if it is not a minhash?
            unimplemented!()
        }
    }
}

#[cfg(test)]
mod test {
    use std::convert::TryInto;
    use std::fs::File;
    use std::io::BufReader;
    use std::path::PathBuf;

    use assert_matches::assert_matches;

    use crate::index::linear::LinearIndex;
    use crate::index::sbt::Factory;
    use crate::index::search::{search_minhashes, search_minhashes_containment};
    use crate::index::storage::ReadData;
    use crate::index::{Index, SigStore, MHMT};
    use crate::signature::Signature;

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn load_mhmt() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/v5_mhmt.sbt.json");

        let sbt = crate::index::MHMT::from_path(filename).expect("Loading error");

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
        let sig_data = sigs[0].clone();

        let leaf: SigStore<_> = sig_data.into();

        let results = sbt.find(search_minhashes, &leaf, 0.5).unwrap();
        //assert_eq!(results.len(), 1);
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

        println!(
            "linear leaves {:?} {:?}",
            linear.datasets.len(),
            linear.datasets
        );

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
    fn roundtrip_mhmt() -> Result<(), Box<dyn std::error::Error>> {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/v5_mhmt.sbt.json");

        let sbt = MHMT::from_path(filename)?;

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

        let mut new_sbt: MHMT = MHMT::builder().storage(None).build();
        let datasets = sbt.signatures();
        for l in datasets {
            new_sbt.insert(l)?;
        }

        /*
        for (i, node) in &sbt.nodes {
            assert_eq!(node.data().unwrap(), new_sbt.nodes[i].data().unwrap());
        }
        */

        assert_eq!(new_sbt.signature_refs().len(), 7);
        println!("new_sbt internal {:?} {:?}", sbt.nodes.len(), sbt.nodes);
        println!("new_sbt leaves {:?} {:?}", sbt.leaves.len(), sbt.leaves);

        let results = new_sbt.find(search_minhashes, &leaf, 0.5)?;
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);
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
}
