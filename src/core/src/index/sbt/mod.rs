pub mod mhbt;

/* FIXME: bring back after boomphf changes
pub mod ukhs;
*/

/* FIXME: bring back after MQF works on macOS and Windows
#[cfg(not(target_arch = "wasm32"))]
pub mod mhmt;
*/

use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::fs::File;
use std::hash::BuildHasherDefault;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use log::info;
use nohash_hasher::NoHashHasher;
use once_cell::sync::OnceCell;
use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

use crate::index::{Comparable, DatasetInfo, Index, SigStore};
use crate::signature::Signature;
use crate::storage::{FSStorage, InnerStorage, ReadData, Storage, StorageInfo};
use crate::Error;

use crate::traits::{ToWriter, Update};

pub trait FromFactory<N> {
    fn factory(&self, name: &str) -> Result<N, Error>;
}

#[derive(TypedBuilder)]
pub struct SBT<N, L> {
    #[builder(default = 2)]
    d: u32,

    #[builder(default, setter(into))]
    storage: Option<InnerStorage>,

    #[builder(default = Factory::GraphFactory { args: (1, 100000.0, 4) })]
    factory: Factory,

    #[builder(default = HashMap::default())]
    nodes: HashMap<u64, N>,

    #[builder(default = HashMap::default())]
    leaves: HashMap<u64, SigStore<L>>,
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

    pub fn storage(&self) -> Option<InnerStorage> {
        self.storage.clone()
    }

    /*
    fn fill_up(&mut self) -> Result<(), Error> {
        let mut visited = HashSet::new();
        let mut queue: Vec<_> = self.leaves.keys().collect();

        while !queue.is_empty() {
            let pos = queue.pop().unwrap();

            if !visited.contains(&pos) {
                visited.insert(pos);
            }
        }

        Ok(())
    }
    */

    // combine
}

impl<T, U> SBT<Node<U>, T>
where
    T: ToWriter + Clone,
    U: ToWriter,
    Node<U>: ReadData<U>,
    SigStore<T>: ReadData<T>,
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

    pub fn from_reader<R, P>(mut rdr: R, path: P) -> Result<SBT<Node<U>, T>, Error>
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
            SBTInfo::V6(ref sbt) => (&sbt.storage.args).into(),
        };
        st.set_base(path.as_ref().to_str().unwrap());
        let storage = InnerStorage::new(st);

        let d = match sinfo {
            SBTInfo::V4(ref sbt) => sbt.d,
            SBTInfo::V5(ref sbt) => sbt.d,
            SBTInfo::V6(ref sbt) => sbt.d,
        };

        let factory = match sinfo {
            SBTInfo::V4(ref sbt) => sbt.factory.clone(),
            SBTInfo::V5(ref sbt) => sbt.factory.clone(),
            SBTInfo::V6(ref sbt) => sbt.factory.clone(),
        };

        let (nodes, leaves) = match sinfo {
            SBTInfo::V6(sbt) => {
                let nodes = sbt
                    .nodes
                    .into_iter()
                    .map(|(n, l)| {
                        (
                            n,
                            Node::builder()
                                .filename(l.filename)
                                .name(l.name)
                                .metadata(l.metadata)
                                .storage(Some(storage.clone()))
                                .build(),
                        )
                    })
                    .collect();
                let leaves = sbt
                    .signatures
                    .into_iter()
                    .map(|(n, l)| {
                        (
                            n,
                            SigStore::builder()
                                .filename(l.filename)
                                .name(l.name)
                                .metadata(l.metadata)
                                .storage(Some(storage.clone()))
                                .build(),
                        )
                    })
                    .collect();
                (nodes, leaves)
            }
            SBTInfo::V5(sbt) => {
                let nodes = sbt
                    .nodes
                    .into_iter()
                    .map(|(n, l)| {
                        (
                            n,
                            Node::builder()
                                .filename(l.filename)
                                .name(l.name)
                                .metadata(l.metadata)
                                .storage(Some(storage.clone()))
                                .build(),
                        )
                    })
                    .collect();
                let leaves = sbt
                    .leaves
                    .into_iter()
                    .map(|(n, l)| {
                        (
                            n,
                            SigStore::builder()
                                .filename(l.filename)
                                .name(l.name)
                                .metadata(l.metadata)
                                .storage(Some(storage.clone()))
                                .build(),
                        )
                    })
                    .collect();
                (nodes, leaves)
            }
            SBTInfo::V4(sbt) => {
                let nodes = sbt
                    .nodes
                    .iter()
                    .filter_map(|(n, x)| match x {
                        NodeInfoV4::Node(l) => Some((
                            *n,
                            Node::builder()
                                .filename(l.filename.clone())
                                .name(l.name.clone())
                                .metadata(l.metadata.clone())
                                .storage(Some(storage.clone()))
                                .build(),
                        )),
                        NodeInfoV4::Leaf(_) => None,
                    })
                    .collect();

                let leaves = sbt
                    .nodes
                    .into_iter()
                    .filter_map(|(n, x)| match x {
                        NodeInfoV4::Node(_) => None,
                        NodeInfoV4::Leaf(l) => Some((
                            n,
                            SigStore::builder()
                                .filename(l.filename)
                                .name(l.name)
                                .metadata(l.metadata)
                                .storage(Some(storage.clone()))
                                .build(),
                        )),
                    })
                    .collect();

                (nodes, leaves)
            }
        };

        Ok(SBT {
            d,
            factory,
            storage: Some(storage),
            nodes,
            leaves,
        })
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<SBT<Node<U>, T>, Error> {
        let file = File::open(&path)?;
        let mut reader = BufReader::new(file);

        // TODO: match with available Storage while we don't
        // add a function to build a Storage from a StorageInfo
        let mut basepath = PathBuf::new();
        basepath.push(path);
        // TODO: canonicalize doesn't work on wasm32-wasi
        //basepath.canonicalize()?;

        let sbt = SBT::<Node<U>, T>::from_reader(&mut reader, &basepath.parent().unwrap())?;
        Ok(sbt)
    }

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
                let subdir = format!(".sbt.{}", basename);
                InnerStorage::new(FSStorage::new(location.to_str().unwrap(), &subdir))
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
                    l.storage = Some(storage.clone());

                    let filename = (*l).save(&l.filename).unwrap();
                    let new_node = NodeInfo {
                        filename,
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
                    l.storage = Some(storage.clone());

                    // TODO: this should be l.md5sum(), not l.filename
                    let filename = (*l).save(&l.filename).unwrap();
                    let new_node = DatasetInfo {
                        filename,
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

    pub fn leaves(&self) -> Vec<SigStore<T>> {
        self.leaves.values().cloned().collect()
    }
}

impl<'a, N, L> Index<'a> for SBT<N, L>
where
    N: Comparable<N> + Comparable<L> + Update<N> + Debug + Default,
    L: Comparable<L> + Update<N> + Clone + Debug + Default,
    SBT<N, L>: FromFactory<N>,
    SigStore<L>: From<L> + ReadData<L>,
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
                    let data = leaf.data().expect("Error reading data");
                    if search_fn(data, sig, threshold) {
                        matches.push(data);
                    }
                }
            }
        }

        Ok(matches)
    }

    fn insert(&mut self, dataset: L) -> Result<(), Error> {
        if self.leaves.is_empty() {
            // in this case the tree is empty,
            // just add the dataset to the first available leaf
            self.leaves.entry(0).or_insert_with(|| dataset.into());
            return Ok(());
        }

        // we can unwrap here because the root node case
        // only happens on an empty tree, and if we got
        // to this point we have at least one leaf already.
        // TODO: find position by similarity search
        let pos = self.leaves.keys().max().unwrap() + 1;
        let parent_pos = self.parent(pos).unwrap();
        let final_pos;

        if let Entry::Occupied(pnode) = self.leaves.entry(parent_pos) {
            // Case 1: parent is a Leaf
            // create a new internal node, add it to self.nodes[parent_pos]

            let (_, leaf) = pnode.remove_entry();

            let mut new_node = self.factory(&format!("internal.{}", parent_pos))?;

            // for each children update the parent node
            // TODO: write the update method
            leaf.data.get().unwrap().update(&mut new_node)?;
            dataset.update(&mut new_node)?;

            // node and parent are children of new internal node
            let mut c_pos = self.children(parent_pos).into_iter().take(2);
            let c1_pos = c_pos.next().unwrap();
            let c2_pos = c_pos.next().unwrap();

            self.leaves.entry(c1_pos).or_insert(leaf);
            self.leaves.entry(c2_pos).or_insert_with(|| dataset.into());
            final_pos = c2_pos;

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
                    self.leaves.entry(pos).or_insert_with(|| dataset.into());
                    final_pos = pos;
                }

                // Case 3: parent is None/empty
                // this can happen with d != 2, need to create parent node
                Entry::Vacant(pnode) => {
                    dataset.update(&mut new_node)?;
                    self.leaves.entry(c_pos).or_insert_with(|| dataset.into());
                    final_pos = c_pos;
                    pnode.insert(new_node);
                }
            }
        }

        let entry = &self.leaves[&final_pos];
        let data = entry.data.get().unwrap();

        let mut parent_pos = parent_pos;
        while let Some(ppos) = self.parent(parent_pos) {
            if let Entry::Occupied(mut pnode) = self.nodes.entry(parent_pos) {
                //TODO: use children for this node to update, instead of dragging
                // dataset up to the root? It would be more generic, but this
                // works for minhash, draff signatures and nodegraphs...
                data.update(&mut pnode.get_mut())?;
            }
            parent_pos = ppos;
        }

        Ok(())
    }

    /*
    fn batch_insert(&mut self, nodes: Vec<Self::Item>) -> Result<(), Error> {
        self = scaffold(nodes, self.storage());
        Ok(())
    }
    */

    fn save<P: AsRef<Path>>(&self, _path: P) -> Result<(), Error> {
        unimplemented!();
    }

    fn load<P: AsRef<Path>>(_path: P) -> Result<(), Error> {
        unimplemented!()
    }

    fn signatures(&self) -> Vec<Self::Item> {
        self.leaves
            .values()
            .map(|x| x.data().unwrap().clone())
            .collect()
    }

    fn signature_refs(&self) -> Vec<&Self::Item> {
        self.leaves.values().map(|x| x.data().unwrap()).collect()
    }

    /*
    fn iter_signatures(&'a self) -> Self::SignatureIterator {
        self.leaves.values()
    }
    */
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
pub struct Node<T> {
    #[builder(setter(into))]
    filename: String,

    #[builder(setter(into))]
    name: String,

    metadata: HashMap<String, u64>,

    #[builder(default)]
    storage: Option<InnerStorage>,

    #[builder(setter(into), default)]
    data: OnceCell<T>,
}

impl<T> Node<T>
where
    T: ToWriter,
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

impl<T> PartialEq for Node<T>
where
    T: PartialEq,
    Node<T>: ReadData<T>,
{
    fn eq(&self, other: &Node<T>) -> bool {
        self.data().unwrap() == other.data().unwrap()
    }
}

impl<T> SigStore<T>
where
    T: ToWriter,
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
    T: Debug,
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
    Leaf(DatasetInfo),
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

#[derive(Serialize, Deserialize)]
struct SBTInfoV6<N, L> {
    d: u32,
    version: u32,
    storage: StorageInfo,
    factory: Factory,
    nodes: HashMap<u64, N>,
    signatures: HashMap<u64, L>,
}

#[derive(Deserialize)]
#[serde(untagged)]
enum SBTInfo {
    V6(SBTInfoV6<NodeInfo, DatasetInfo>),
    V5(SBTInfoV5<NodeInfo, DatasetInfo>),
    V4(SBTInfoV4<NodeInfoV4>),
}

enum BinaryTree {
    Empty,
    Internal(Box<TreeNode<HashSet<u64, BuildHasherDefault<NoHashHasher<u64>>>>>),
    Leaf(Box<TreeNode<SigStore<Signature>>>),
}

struct TreeNode<T> {
    element: T,
    left: BinaryTree,
    right: BinaryTree,
}

pub fn scaffold<N>(
    mut datasets: Vec<SigStore<Signature>>,
    storage: Option<InnerStorage>,
) -> SBT<Node<N>, Signature>
where
    N: Clone + Default,
{
    let mut leaves: HashMap<u64, SigStore<Signature>> = HashMap::with_capacity(datasets.len());

    let mut next_round = Vec::new();

    // generate two bottom levels:
    // - datasets
    // - first level of internal nodes
    info!("Start processing leaves");
    while !datasets.is_empty() {
        let next_leaf = datasets.pop().unwrap();

        let (simleaf_tree, in_common) = if datasets.is_empty() {
            (BinaryTree::Empty, next_leaf.mins().into_iter().collect())
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

            let in_common = next_leaf
                .mins()
                .into_iter()
                .collect::<HashSet<u64, BuildHasherDefault<NoHashHasher<u64>>>>()
                .union(&similar_leaf.mins().into_iter().collect())
                .cloned()
                .collect();

            let simleaf_tree = BinaryTree::Leaf(Box::new(TreeNode {
                element: similar_leaf,
                left: BinaryTree::Empty,
                right: BinaryTree::Empty,
            }));
            (simleaf_tree, in_common)
        };

        let leaf_tree = BinaryTree::Leaf(Box::new(TreeNode {
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
            info!("Processed {} leaves", next_round.len() * 2);
        }
    }
    info!("Finished processing leaves");

    // while we don't get to the root, generate intermediary levels
    while next_round.len() != 1 {
        next_round = BinaryTree::process_internal_level(next_round);
        info!("Finished processing round {}", next_round.len());
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
                BinaryTree::Leaf(leaf) => {
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

    // Remove this when MSRV is >= 1.40
    #[allow(clippy::mem_replace_with_default)]
    fn new_tree(mut left: BinaryTree, mut right: BinaryTree) -> BinaryTree {
        let in_common = if let BinaryTree::Internal(ref mut el1) = left {
            match right {
                BinaryTree::Internal(ref mut el2) => {
                    let c1 = std::mem::replace(
                        &mut el1.element,
                        HashSet::<u64, BuildHasherDefault<NoHashHasher<u64>>>::default(),
                    );
                    let c2 = std::mem::replace(
                        &mut el2.element,
                        HashSet::<u64, BuildHasherDefault<NoHashHasher<u64>>>::default(),
                    );
                    c1.union(&c2).cloned().collect()
                }
                BinaryTree::Empty => std::mem::replace(
                    &mut el1.element,
                    HashSet::<u64, BuildHasherDefault<NoHashHasher<u64>>>::default(),
                ),
                _ => panic!("Should not see a Leaf at this level"),
            }
        } else {
            HashSet::<u64, BuildHasherDefault<NoHashHasher<u64>>>::default()
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

/*
impl<U> From<LinearIndex<Signature>> for SBT<Node<U>, Signature>
where
    U: Default + Clone,
{
    fn from(other: LinearIndex<Signature>) -> Self {
        let storage = other.storage();
        scaffold(other.datasets, storage)
    }
}
*/
