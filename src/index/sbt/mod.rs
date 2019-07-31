pub mod mhbt;
pub mod ukhs;

use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::fs::File;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{BufReader, Read};
use std::iter::FromIterator;
use std::mem;
use std::path::{Path, PathBuf};
use std::rc::Rc;

use failure::Error;
use lazy_init::Lazy;
use serde_derive::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

use crate::index::storage::{FSStorage, ReadData, Storage, StorageInfo, ToWriter};
use crate::index::{Comparable, Dataset, DatasetInfo, Index};
use crate::signature::Signature;

pub trait Update<O> {
    fn update(&self, other: &mut O) -> Result<(), Error>;
}

pub trait FromFactory<N> {
    fn factory(&self, name: &str) -> Result<N, Error>;
}

#[derive(TypedBuilder)]
pub struct SBT<N, L> {
    #[builder(default = 2)]
    d: u32,

    #[builder(default)]
    storage: Option<Rc<dyn Storage>>,

    #[builder(default_code = r#"Factory::GraphFactory { args: (1, 100000.0, 4) }"#)]
    factory: Factory,

    #[builder(default_code = "HashMap::default()")]
    nodes: HashMap<u64, N>,

    #[builder(default_code = "HashMap::default()")]
    leaves: HashMap<u64, L>,
}

const fn parent(pos: u64, d: u64) -> u64 {
    ((pos - 1) / d) as u64
}

const fn child(parent: u64, pos: u64, d: u64) -> u64 {
    d * parent + pos + 1
}

impl<N, L> SBT<N, L>
where
    L: std::clone::Clone + Default,
    N: Default,
{
    #[inline(always)]
    fn parent(&self, pos: u64) -> Option<u64> {
        if pos == 0 {
            None
        } else {
            Some(parent(pos, u64::from(self.d)))
        }
    }

    #[inline(always)]
    fn child(&self, parent: u64, pos: u64) -> u64 {
        child(parent, pos, u64::from(self.d))
    }

    #[inline(always)]
    fn children(&self, pos: u64) -> Vec<u64> {
        (0..u64::from(self.d)).map(|c| self.child(pos, c)).collect()
    }

    pub fn storage(&self) -> Option<Rc<dyn Storage>> {
        self.storage.clone()
    }

    // combine
}

impl<T, U> SBT<Node<U>, Dataset<T>>
where
    T: std::marker::Sync + ToWriter,
    U: std::marker::Sync + ToWriter,
    Node<U>: ReadData<U>,
    Dataset<T>: ReadData<T>,
{
    fn parse_v4<R>(rdr: &mut R) -> Result<SBTInfo, Error>
    where
        R: Read,
    {
        let sinfo: SBTInfoV4<NodeInfoV4> = serde_json::from_reader(rdr)?;
        Ok(SBTInfo::V4(sinfo))
    }

    fn parse_v5<R>(rdr: &mut R) -> Result<SBTInfo, Error>
    where
        R: Read,
    {
        let sinfo: SBTInfoV5<NodeInfo, DatasetInfo> = serde_json::from_reader(rdr)?;
        Ok(SBTInfo::V5(sinfo))
    }

    pub fn from_reader<R, P>(rdr: &mut R, path: P) -> Result<SBT<Node<U>, Dataset<T>>, Error>
    where
        R: Read,
        P: AsRef<Path>,
    {
        // TODO: I would love to do this, but I get an untagged enum error with
        // SBTInfo...
        //let sinfo: SBTInfo = serde_json::from_reader(rdr)?;

        let mut s = String::new();
        rdr.read_to_string(&mut s)?;

        let sinfo =
            Self::parse_v5(&mut s.as_bytes()).or_else(|_| Self::parse_v4(&mut s.as_bytes()))?;

        // TODO: support other storages
        let mut st: FSStorage = match sinfo {
            SBTInfo::V4(ref sbt) => (&sbt.storage.args).into(),
            SBTInfo::V5(ref sbt) => (&sbt.storage.args).into(),
        };
        st.set_base(path.as_ref().to_str().unwrap());
        let storage: Rc<dyn Storage> = Rc::new(st);

        let d = match sinfo {
            SBTInfo::V4(ref sbt) => sbt.d,
            SBTInfo::V5(ref sbt) => sbt.d,
        };

        let factory = match sinfo {
            SBTInfo::V4(ref sbt) => sbt.factory.clone(),
            SBTInfo::V5(ref sbt) => sbt.factory.clone(),
        };

        let (nodes, leaves) = match sinfo {
            SBTInfo::V5(sbt) => {
                let nodes = sbt
                    .nodes
                    .into_iter()
                    .map(|(n, l)| {
                        let new_node = Node {
                            filename: l.filename,
                            name: l.name,
                            metadata: l.metadata,
                            storage: Some(Rc::clone(&storage)),
                            data: Rc::new(Lazy::new()),
                        };
                        (n, new_node)
                    })
                    .collect();
                let leaves = sbt
                    .leaves
                    .into_iter()
                    .map(|(n, l)| {
                        let new_node = Dataset {
                            filename: l.filename,
                            name: l.name,
                            metadata: l.metadata,
                            storage: Some(Rc::clone(&storage)),
                            data: Rc::new(Lazy::new()),
                        };
                        (n, new_node)
                    })
                    .collect();
                (nodes, leaves)
            }
            SBTInfo::V4(sbt) => {
                let nodes = sbt
                    .nodes
                    .iter()
                    .filter_map(|(n, x)| match x {
                        NodeInfoV4::Node(l) => {
                            let new_node = Node {
                                filename: l.filename.clone(),
                                name: l.name.clone(),
                                metadata: l.metadata.clone(),
                                storage: Some(Rc::clone(&storage)),
                                data: Rc::new(Lazy::new()),
                            };
                            Some((*n, new_node))
                        }
                        NodeInfoV4::Dataset(_) => None,
                    })
                    .collect();

                let leaves = sbt
                    .nodes
                    .into_iter()
                    .filter_map(|(n, x)| match x {
                        NodeInfoV4::Node(_) => None,
                        NodeInfoV4::Dataset(l) => {
                            let new_node = Dataset {
                                filename: l.filename,
                                name: l.name,
                                metadata: l.metadata,
                                storage: Some(Rc::clone(&storage)),
                                data: Rc::new(Lazy::new()),
                            };
                            Some((n, new_node))
                        }
                    })
                    .collect();

                (nodes, leaves)
            }
        };

        Ok(SBT {
            d,
            factory,
            storage: Some(Rc::clone(&storage)),
            nodes,
            leaves,
        })
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<SBT<Node<U>, Dataset<T>>, Error> {
        let file = File::open(&path)?;
        let mut reader = BufReader::new(file);

        // TODO: match with available Storage while we don't
        // add a function to build a Storage from a StorageInfo
        let mut basepath = PathBuf::new();
        basepath.push(path);
        basepath.canonicalize()?;

        let sbt =
            SBT::<Node<U>, Dataset<T>>::from_reader(&mut reader, &basepath.parent().unwrap())?;
        Ok(sbt)
    }

    pub fn save_file<P: AsRef<Path>>(
        &mut self,
        path: P,
        storage: Option<Rc<dyn Storage>>,
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
                let subdir = format!(".sbt.{}", basename);
                Rc::new(FSStorage::new(location.to_str().unwrap(), &subdir))
            }
        };

        let args = storage.args();
        let storage_info = StorageInfo {
            backend: "FSStorage".into(),
            args,
        };

        let info: SBTInfoV5<NodeInfo, DatasetInfo> = SBTInfoV5 {
            d: self.d,
            factory: self.factory.clone(),
            storage: storage_info,
            version: 5,
            nodes: self
                .nodes
                .iter_mut()
                .map(|(n, l)| {
                    // Trigger data loading
                    let _: &U = (*l).data().expect("Couldn't load data");

                    // set storage to new one
                    mem::replace(&mut l.storage, Some(Rc::clone(&storage)));

                    let filename = (*l).save(&l.filename).unwrap();
                    let new_node = NodeInfo {
                        filename: filename,
                        name: l.name.clone(),
                        metadata: l.metadata.clone(),
                    };
                    (*n, new_node)
                })
                .collect(),
            leaves: self
                .leaves
                .iter_mut()
                .map(|(n, l)| {
                    // Trigger data loading
                    let _: &T = (*l).data().unwrap();

                    // set storage to new one
                    mem::replace(&mut l.storage, Some(Rc::clone(&storage)));

                    let filename = (*l).save(&l.filename).unwrap();
                    let new_node = DatasetInfo {
                        filename: filename,
                        name: l.name.clone(),
                        metadata: l.metadata.clone(),
                    };
                    (*n, new_node)
                })
                .collect(),
        };

        let file = File::create(path)?;
        serde_json::to_writer(file, &info)?;

        Ok(())
    }
}

impl<N, L> Index for SBT<N, L>
where
    N: Comparable<N> + Comparable<L> + Update<N> + Debug + Default,
    L: Comparable<L> + Update<N> + Clone + Debug + Default,
    SBT<N, L>: FromFactory<N>,
{
    type Item = L;

    fn find<F>(&self, search_fn: F, sig: &L, threshold: f64) -> Result<Vec<&L>, Error>
    where
        F: Fn(&dyn Comparable<Self::Item>, &Self::Item, f64) -> bool,
    {
        let mut matches = Vec::new();
        let mut visited = HashSet::new();
        let mut queue = vec![0u64];

        while !queue.is_empty() {
            let pos = queue.pop().unwrap();
            if !visited.contains(&pos) {
                visited.insert(pos);

                if let Some(node) = self.nodes.get(&pos) {
                    if search_fn(&node, sig, threshold) {
                        for c in self.children(pos) {
                            queue.push(c);
                        }
                    }
                } else if let Some(leaf) = self.leaves.get(&pos) {
                    if search_fn(leaf, sig, threshold) {
                        matches.push(leaf);
                    }
                }
            }
        }

        Ok(matches)
    }

    fn insert(&mut self, dataset: &L) -> Result<(), Error> {
        if self.leaves.is_empty() {
            // in this case the tree is empty,
            // just add the dataset to the first available leaf
            self.leaves.entry(0).or_insert(dataset.clone());
            return Ok(());
        }

        // we can unwrap here because the root node case
        // only happens on an empty tree, and if we got
        // to this point we have at least one leaf already.
        // TODO: find position by similarity search
        let pos = self.leaves.keys().max().unwrap() + 1;
        let parent_pos = self.parent(pos).unwrap();

        if let Entry::Occupied(pnode) = self.leaves.entry(parent_pos) {
            // Case 1: parent is a Leaf
            // create a new internal node, add it to self.nodes[parent_pos]

            let (_, leaf) = pnode.remove_entry();

            let mut new_node = self.factory(&format!("internal.{}", parent_pos))?;

            // for each children update the parent node
            // TODO: write the update method
            leaf.update(&mut new_node)?;
            dataset.update(&mut new_node)?;

            // node and parent are children of new internal node
            let mut c_pos = self.children(parent_pos).into_iter().take(2);
            let c1_pos = c_pos.next().unwrap();
            let c2_pos = c_pos.next().unwrap();

            self.leaves.entry(c1_pos).or_insert(leaf);
            self.leaves.entry(c2_pos).or_insert(dataset.clone());

            // add the new internal node to self.nodes[parent_pos)
            // TODO check if it is really empty?
            self.nodes.entry(parent_pos).or_insert(new_node);
        } else {
            // TODO: moved these two lines here to avoid borrow checker
            // error E0502 in the Vacant case, but would love to avoid it!
            let mut new_node = self.factory(&format!("internal.{}", parent_pos))?;
            let c_pos = self.children(parent_pos)[0];

            match self.nodes.entry(parent_pos) {
                // Case 2: parent is a node and has an empty child spot available
                // (if there isn't an empty spot, it was already covered by case 1)
                Entry::Occupied(mut pnode) => {
                    dataset.update(&mut pnode.get_mut())?;
                    self.leaves.entry(pos).or_insert(dataset.clone());
                }

                // Case 3: parent is None/empty
                // this can happen with d != 2, need to create parent node
                Entry::Vacant(pnode) => {
                    self.leaves.entry(c_pos).or_insert(dataset.clone());
                    dataset.update(&mut new_node)?;
                    pnode.insert(new_node);
                }
            }
        }

        let mut parent_pos = parent_pos;
        while let Some(ppos) = self.parent(parent_pos) {
            if let Entry::Occupied(mut pnode) = self.nodes.entry(parent_pos) {
                //TODO: use children for this node to update, instead of dragging
                // dataset up to the root? It would be more generic, but this
                // works for minhash, draff signatures and nodegraphs...
                dataset.update(&mut pnode.get_mut())?;
            }
            parent_pos = ppos;
        }

        Ok(())
    }

    fn save<P: AsRef<Path>>(&self, _path: P) -> Result<(), Error> {
        unimplemented!()
    }

    fn load<P: AsRef<Path>>(_path: P) -> Result<(), Error> {
        unimplemented!()
    }

    fn datasets(&self) -> Vec<Self::Item> {
        self.leaves.values().cloned().collect()
    }
}

/*
#[derive(TypedBuilder, Clone, Default, Serialize, Deserialize)]
pub struct Factory {
    class: String,
    args: Vec<u64>,
}
*/

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "class")]
pub enum Factory {
    GraphFactory { args: (u64, f64, u64) },
}

#[derive(TypedBuilder, Default, Clone)]
pub struct Node<T>
where
    T: std::marker::Sync,
{
    filename: String,
    name: String,
    metadata: HashMap<String, u64>,
    storage: Option<Rc<dyn Storage>>,
    #[builder(default)]
    pub(crate) data: Rc<Lazy<T>>,
}

impl<T> Node<T>
where
    T: Sync + ToWriter,
{
    pub fn save(&self, path: &str) -> Result<String, Error> {
        if let Some(storage) = &self.storage {
            if let Some(data) = self.data.get() {
                let mut buffer = Vec::new();
                data.to_writer(&mut buffer)?;

                Ok(storage.save(path, &buffer)?)
            } else {
                // TODO throw error, data was not initialized
                unimplemented!()
            }
        } else {
            unimplemented!()
        }
    }
}

impl<T> Dataset<T>
where
    T: Sync + ToWriter,
{
    pub fn save(&self, path: &str) -> Result<String, Error> {
        if let Some(storage) = &self.storage {
            if let Some(data) = self.data.get() {
                let mut buffer = Vec::new();
                data.to_writer(&mut buffer)?;

                Ok(storage.save(path, &buffer)?)
            } else {
                unimplemented!()
            }
        } else {
            unimplemented!()
        }
    }
}

impl<T> std::fmt::Debug for Node<T>
where
    T: std::marker::Sync + std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Node [name={}, filename={}, metadata: {:?}, data: {:?}]",
            self.name,
            self.filename,
            self.metadata,
            self.data.get().is_some()
        )
    }
}

#[derive(Serialize, Deserialize, Debug)]
struct NodeInfo {
    filename: String,
    name: String,
    metadata: HashMap<String, u64>,
}

#[derive(Serialize, Deserialize, Debug)]
#[serde(untagged)]
enum NodeInfoV4 {
    Node(NodeInfo),
    Dataset(DatasetInfo),
}

#[derive(Serialize, Deserialize)]
struct SBTInfoV4<N> {
    d: u32,
    version: u32,
    storage: StorageInfo,
    factory: Factory,
    nodes: HashMap<u64, N>,
}

#[derive(Serialize, Deserialize)]
struct SBTInfoV5<N, L> {
    d: u32,
    version: u32,
    storage: StorageInfo,
    factory: Factory,
    nodes: HashMap<u64, N>,
    leaves: HashMap<u64, L>,
}

#[derive(Deserialize)]
#[serde(untagged)]
enum SBTInfo {
    V5(SBTInfoV5<NodeInfo, DatasetInfo>),
    V4(SBTInfoV4<NodeInfoV4>),
}

// This comes from finch
pub struct NoHashHasher(u64);

impl Default for NoHashHasher {
    #[inline]
    fn default() -> NoHashHasher {
        NoHashHasher(0x0)
    }
}

impl Hasher for NoHashHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        *self = NoHashHasher(
            (u64::from(bytes[0]) << 24)
                + (u64::from(bytes[1]) << 16)
                + (u64::from(bytes[2]) << 8)
                + u64::from(bytes[3]),
        );
    }
    fn finish(&self) -> u64 {
        self.0
    }
}

type HashIntersection = HashSet<u64, BuildHasherDefault<NoHashHasher>>;

enum BinaryTree {
    Empty,
    Internal(Box<TreeNode<HashIntersection>>),
    Dataset(Box<TreeNode<Dataset<Signature>>>),
}

struct TreeNode<T> {
    element: T,
    left: BinaryTree,
    right: BinaryTree,
}

pub fn scaffold<N>(
    mut datasets: Vec<Dataset<Signature>>,
    storage: Option<Rc<dyn Storage>>,
) -> SBT<Node<N>, Dataset<Signature>>
where
    N: std::marker::Sync + std::clone::Clone + std::default::Default,
{
    let mut leaves: HashMap<u64, Dataset<Signature>> = HashMap::with_capacity(datasets.len());

    let mut next_round = Vec::new();

    // generate two bottom levels:
    // - datasets
    // - first level of internal nodes
    eprintln!("Start processing leaves");
    while !datasets.is_empty() {
        let next_leaf = datasets.pop().unwrap();

        let (simleaf_tree, in_common) = if datasets.is_empty() {
            (
                BinaryTree::Empty,
                HashIntersection::from_iter(next_leaf.mins().into_iter()),
            )
        } else {
            let mut similar_leaf_pos = 0;
            let mut current_max = 0;
            for (pos, leaf) in datasets.iter().enumerate() {
                let common = next_leaf.count_common(leaf);
                if common > current_max {
                    current_max = common;
                    similar_leaf_pos = pos;
                }
            }

            let similar_leaf = datasets.remove(similar_leaf_pos);

            let in_common = HashIntersection::from_iter(next_leaf.mins().into_iter())
                .union(&HashIntersection::from_iter(
                    similar_leaf.mins().into_iter(),
                ))
                .cloned()
                .collect();

            let simleaf_tree = BinaryTree::Dataset(Box::new(TreeNode {
                element: similar_leaf,
                left: BinaryTree::Empty,
                right: BinaryTree::Empty,
            }));
            (simleaf_tree, in_common)
        };

        let leaf_tree = BinaryTree::Dataset(Box::new(TreeNode {
            element: next_leaf,
            left: BinaryTree::Empty,
            right: BinaryTree::Empty,
        }));

        let tree = BinaryTree::Internal(Box::new(TreeNode {
            element: in_common,
            left: leaf_tree,
            right: simleaf_tree,
        }));

        next_round.push(tree);

        if next_round.len() % 100 == 0 {
            eprintln!("Processed {} leaves", next_round.len() * 2);
        }
    }
    eprintln!("Finished processing leaves");

    // while we don't get to the root, generate intermediary levels
    while next_round.len() != 1 {
        next_round = BinaryTree::process_internal_level(next_round);
        eprintln!("Finished processing round {}", next_round.len());
    }

    // Convert from binary tree to nodes/leaves
    let root = next_round.pop().unwrap();
    let mut visited = HashSet::new();
    let mut queue = vec![(0u64, root)];

    while !queue.is_empty() {
        let (pos, cnode) = queue.pop().unwrap();
        if !visited.contains(&pos) {
            visited.insert(pos);

            match cnode {
                BinaryTree::Dataset(leaf) => {
                    leaves.insert(pos, leaf.element);
                }
                BinaryTree::Internal(mut node) => {
                    let left = std::mem::replace(&mut node.left, BinaryTree::Empty);
                    let right = std::mem::replace(&mut node.right, BinaryTree::Empty);
                    queue.push((2 * pos + 1, left));
                    queue.push((2 * pos + 2, right));
                }
                BinaryTree::Empty => (),
            }
        }
    }

    SBT::builder()
        .storage(storage)
        .nodes(HashMap::default())
        .leaves(leaves)
        .build()
}

impl BinaryTree {
    fn process_internal_level(mut current_round: Vec<BinaryTree>) -> Vec<BinaryTree> {
        let mut next_round = Vec::with_capacity(current_round.len() + 1);

        while !current_round.is_empty() {
            let next_node = current_round.pop().unwrap();

            let similar_node = if current_round.is_empty() {
                BinaryTree::Empty
            } else {
                let mut similar_node_pos = 0;
                let mut current_max = 0;
                for (pos, cmpe) in current_round.iter().enumerate() {
                    let common = BinaryTree::intersection_size(&next_node, &cmpe);
                    if common > current_max {
                        current_max = common;
                        similar_node_pos = pos;
                    }
                }
                current_round.remove(similar_node_pos)
            };

            let tree = BinaryTree::new_tree(next_node, similar_node);

            next_round.push(tree);
        }
        next_round
    }

    fn new_tree(mut left: BinaryTree, mut right: BinaryTree) -> BinaryTree {
        let in_common = if let BinaryTree::Internal(ref mut el1) = left {
            match right {
                BinaryTree::Internal(ref mut el2) => {
                    let c1 = std::mem::replace(&mut el1.element, HashIntersection::default());
                    let c2 = std::mem::replace(&mut el2.element, HashIntersection::default());
                    c1.union(&c2).cloned().collect()
                }
                BinaryTree::Empty => {
                    std::mem::replace(&mut el1.element, HashIntersection::default())
                }
                _ => panic!("Should not see a Dataset at this level"),
            }
        } else {
            HashIntersection::default()
        };

        BinaryTree::Internal(Box::new(TreeNode {
            element: in_common,
            left,
            right,
        }))
    }

    fn intersection_size(n1: &BinaryTree, n2: &BinaryTree) -> usize {
        if let BinaryTree::Internal(ref el1) = n1 {
            if let BinaryTree::Internal(ref el2) = n2 {
                return el1.element.intersection(&el2.element).count();
            }
        };
        0
    }
}

#[cfg(test)]
mod test {
    use std::fs::File;
    use std::io::{BufReader, Seek, SeekFrom};
    use std::path::PathBuf;
    use std::rc::Rc;
    use tempfile;

    use assert_matches::assert_matches;
    use lazy_init::Lazy;

    use super::{scaffold, Factory};

    use crate::index::linear::LinearIndex;
    use crate::index::search::{search_minhashes, search_minhashes_containment};
    use crate::index::{Dataset, Index, MHBT};
    use crate::signature::Signature;

    #[test]
    fn save_sbt() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/test-data/v5.sbt.json");

        let mut sbt = MHBT::from_path(filename).expect("Loading error");

        let mut tmpfile = tempfile::NamedTempFile::new().unwrap();
        sbt.save_file(tmpfile.path(), None).unwrap();

        tmpfile.seek(SeekFrom::Start(0)).unwrap();
    }

    #[test]
    fn load_sbt() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/test-data/v5.sbt.json");

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
        filename.push("tests/test-data/.sbt.v3/60f7e23c24a8d94791cc7a8680c493f9");

        let mut reader = BufReader::new(File::open(filename).unwrap());
        let sigs = Signature::load_signatures(&mut reader, 31, Some("DNA".into()), None).unwrap();
        let sig_data = sigs[0].clone();

        let data = Lazy::new();
        data.get_or_create(|| sig_data);

        let leaf = Dataset::builder()
            .data(Rc::new(data))
            .filename("")
            .name("")
            .metadata("")
            .storage(None)
            .build();

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
            linear.insert(l.1).unwrap();
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
    fn scaffold_sbt() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/test-data/v5.sbt.json");

        let sbt = MHBT::from_path(filename).expect("Loading error");

        let new_sbt: MHBT = scaffold(sbt.datasets(), sbt.storage());

        assert_eq!(new_sbt.datasets().len(), 7);
    }

    #[test]
    fn load_v4() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/test-data/v4.sbt.json");

        let _sbt = MHBT::from_path(filename).expect("Loading error");
    }

    #[test]
    fn load_v5() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/test-data/v5.sbt.json");

        let _sbt = MHBT::from_path(filename).expect("Loading error");
    }
}
