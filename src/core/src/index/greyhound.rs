use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use cfg_if::cfg_if;
use getset::{CopyGetters, Getters, Setters};
use log::{debug, info};
use nohash_hasher::BuildNoHashHasher;
use serde::{Deserialize, Serialize};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::encodings::{Color, Colors, Idx};
use crate::index::Index;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::minhash::KmerMinHash;
use crate::sketch::Sketch;
//use crate::storage::{InnerStorage, MemStorage};
use crate::Error;
use crate::HashIntoType;

type HashToColor = HashMap<HashIntoType, Color, BuildNoHashHasher<HashIntoType>>;
type SigCounter = counter::Counter<Idx>;

#[derive(Serialize, Deserialize)]
pub struct RevIndex {
    hash_to_color: HashToColor,

    sig_files: Vec<PathBuf>,

    #[serde(skip)]
    ref_sigs: Option<Vec<Signature>>,

    template: Sketch,
    colors: Colors,
    //#[serde(skip)]
    //storage: Option<InnerStorage>,
}

impl RevIndex {
    pub fn load<P: AsRef<Path>>(
        index_path: P,
        queries: Option<&[KmerMinHash]>,
    ) -> Result<RevIndex, Box<dyn std::error::Error>> {
        let (rdr, _) = niffler::from_path(index_path)?;
        let revindex = if let Some(qs) = queries {
            // TODO: avoid loading full revindex if query != None
            /*
            struct PartialRevIndex<T> {
                hashes_to_keep: Option<HashSet<HashIntoType>>,
                marker: PhantomData<fn() -> T>,
            }

            impl<T> PartialRevIndex<T> {
                pub fn new(hashes_to_keep: HashSet<u64>) -> Self {
                    PartialRevIndex {
                        hashes_to_keep: Some(hashes_to_keep),
                        marker: PhantomData,
                    }
                }
            }
            */

            let mut hashes: HashSet<u64> = HashSet::new();
            for q in qs {
                hashes.extend(q.iter_mins());
            }

            //let mut revindex: RevIndex = PartialRevIndex::new(hashes).deserialize(&rdr).unwrap();

            let mut revindex: RevIndex = serde_json::from_reader(rdr)?;
            revindex
                .hash_to_color
                .retain(|hash, _| hashes.contains(hash));
            revindex
        } else {
            // Load the full revindex
            serde_json::from_reader(rdr)?
        };

        Ok(revindex)
    }

    pub fn new(
        search_sigs: &[PathBuf],
        template: &Sketch,
        threshold: usize,
        queries: Option<&[KmerMinHash]>,
        keep_sigs: bool,
    ) -> RevIndex {
        // If threshold is zero, let's merge all queries and save time later
        let merged_query = if let Some(qs) = queries {
            if threshold == 0 {
                let mut merged = qs[0].clone();
                for query in &qs[1..] {
                    merged.merge(query).unwrap();
                }
                Some(merged)
            } else {
                None
            }
        } else {
            None
        };

        let processed_sigs = AtomicUsize::new(0);

        cfg_if! {
          if #[cfg(feature = "parallel")] {
            let (hash_to_color, colors) = search_sigs
                .par_iter()
                .enumerate()
                .filter_map(|(dataset_id, filename)| {
                    let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
                    if i % 1000 == 0 {
                        info!("Processed {} reference sigs", i);
                    }

                    RevIndex::map_hashes_colors(
                        dataset_id,
                        filename,
                        queries,
                        &merged_query,
                        threshold,
                        template,
                    )
                })
                .reduce(
                    || {
                        (
                            HashToColor::with_hasher(BuildNoHashHasher::default()),
                            Colors::default(),
                        )
                    },
                    RevIndex::reduce_hashes_colors,
                );
          } else {
            let (hash_to_color, colors) = search_sigs
                .iter()
                .enumerate()
                .filter_map(|(dataset_id, filename)| {
                    let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
                    if i % 1000 == 0 {
                        info!("Processed {} reference sigs", i);
                    }
                    RevIndex::map_hashes_colors(
                        dataset_id,
                        filename,
                        queries,
                        &merged_query,
                        threshold,
                        template,
                    )
                })
                .fold(
                    (
                        HashToColor::with_hasher(BuildNoHashHasher::default()),
                        Colors::default(),
                    ),
                    RevIndex::reduce_hashes_colors,
                );
          }
        }

        // TODO: build this together with hash_to_idx?
        // TODO: do a parallel (par_iter) and sequential (iter) version
        let ref_sigs = if keep_sigs {
            cfg_if! {
              if #[cfg(feature = "parallel")] {
                Some(
                    search_sigs
                        .par_iter()
                        .map(|ref_path| {
                            Signature::from_path(&ref_path)
                                .unwrap_or_else(|_| panic!("Error processing {:?}", ref_path))
                                .swap_remove(0)
                        })
                        .collect(),
                )
              } else {
                Some(
                    search_sigs
                        .iter()
                        .map(|ref_path| {
                            Signature::from_path(&ref_path)
                                .unwrap_or_else(|_| panic!("Error processing {:?}", ref_path))
                                .swap_remove(0)
                        })
                        .collect(),
                )
              }
            }
        } else {
            None
        };

        RevIndex {
            hash_to_color,
            sig_files: search_sigs.into(),
            ref_sigs,
            template: template.clone(),
            colors,
            //            storage: Some(InnerStorage::new(MemStorage::default())),
        }
    }

    pub fn new_with_sigs(
        search_sigs: Vec<Signature>,
        template: &Sketch,
        threshold: usize,
        queries: Option<&[KmerMinHash]>,
    ) -> RevIndex {
        // If threshold is zero, let's merge all queries and save time later
        let merged_query = if let Some(qs) = queries {
            if threshold == 0 {
                let mut merged = qs[0].clone();
                for query in &qs[1..] {
                    merged.merge(query).unwrap();
                }
                Some(merged)
            } else {
                None
            }
        } else {
            None
        };

        let processed_sigs = AtomicUsize::new(0);

        cfg_if! {
          if #[cfg(feature = "parallel")] {
            let (hash_to_color, colors) = search_sigs
                .par_iter()
                .enumerate()
                .filter_map(|(dataset_id, sig)| {
                    let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
                    if i % 1000 == 0 {
                        info!("Processed {} reference sigs", i);
                    }

                    RevIndex::map_hashes_colors_sigs(
                        dataset_id,
                        sig,
                        queries,
                        &merged_query,
                        threshold,
                        template,
                    )
                })
                .reduce(
                    || {
                        (
                            HashToColor::with_hasher(BuildNoHashHasher::default()),
                            Colors::default(),
                        )
                    },
                    RevIndex::reduce_hashes_colors,
                );
          } else {
            let (hash_to_color, colors) = search_sigs
                .iter()
                .enumerate()
                .filter_map(|(dataset_id, filename)| {
                    let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
                    if i % 1000 == 0 {
                        info!("Processed {} reference sigs", i);
                    }
                    RevIndex::map_hashes_colors_sigs(
                        dataset_id,
                        filename,
                        queries,
                        &merged_query,
                        threshold,
                        template,
                    )
                })
                .fold(
                    (
                        HashToColor::with_hasher(BuildNoHashHasher::default()),
                        Colors::default(),
                    ),
                    RevIndex::reduce_hashes_colors,
                );
          }
        }

        RevIndex {
            hash_to_color,
            sig_files: vec![],
            ref_sigs: search_sigs.into(),
            template: template.clone(),
            colors,
            //storage: None,
        }
    }

    fn map_hashes_colors_sigs(
        dataset_id: usize,
        search_sig: &Signature,
        queries: Option<&[KmerMinHash]>,
        merged_query: &Option<KmerMinHash>,
        threshold: usize,
        template: &Sketch,
    ) -> Option<(HashToColor, Colors)> {
        let mut search_mh = None;
        if let Some(Sketch::MinHash(mh)) = search_sig.select_sketch(&template) {
            search_mh = Some(mh);
        }

        let search_mh = search_mh.expect("Couldn't find a compatible MinHash");
        let mut hash_to_color = HashToColor::with_hasher(BuildNoHashHasher::default());
        let mut colors = Colors::default();
        let mut color = None;

        let mut add_to = |matched_hashes: Vec<u64>, intersection| {
            if !matched_hashes.is_empty() || intersection > threshold as u64 {
                matched_hashes.into_iter().for_each(|hash| {
                    color = Some(colors.update(color, &[dataset_id as Idx]).unwrap());
                    hash_to_color.insert(hash, color.unwrap());
                });
            }
        };

        if let Some(qs) = queries {
            if let Some(ref merged) = merged_query {
                let (matched_hashes, intersection) = merged.intersection(search_mh).unwrap();
                add_to(matched_hashes, intersection);
            } else {
                for query in qs {
                    let (matched_hashes, intersection) = query.intersection(search_mh).unwrap();
                    add_to(matched_hashes, intersection);
                }
            }
        } else {
            let matched = search_mh.mins();
            let size = matched.len() as u64;
            add_to(matched, size);
        };

        if hash_to_color.is_empty() {
            None
        } else {
            Some((hash_to_color, colors))
        }
    }

    fn map_hashes_colors(
        dataset_id: usize,
        filename: &Path,
        queries: Option<&[KmerMinHash]>,
        merged_query: &Option<KmerMinHash>,
        threshold: usize,
        template: &Sketch,
    ) -> Option<(HashToColor, Colors)> {
        let search_sig = Signature::from_path(&filename)
            .unwrap_or_else(|_| panic!("Error processing {:?}", filename))
            .swap_remove(0);

        RevIndex::map_hashes_colors_sigs(
            dataset_id,
            &search_sig,
            queries,
            merged_query,
            threshold,
            template,
        )
    }

    fn reduce_hashes_colors(
        a: (HashToColor, Colors),
        b: (HashToColor, Colors),
    ) -> (HashToColor, Colors) {
        let ((small_hashes, small_colors), (mut large_hashes, mut large_colors)) =
            if a.0.len() > b.0.len() {
                (b, a)
            } else {
                (a, b)
            };

        small_hashes.into_iter().for_each(|(hash, color)| {
            large_hashes
                .entry(hash)
                .and_modify(|entry| {
                    // Hash is already present.
                    // Update the current color by adding the indices from
                    // small_colors.
                    let ids = small_colors.indices(&color);
                    let new_color = large_colors.update(Some(*entry), ids).unwrap();
                    *entry = new_color;
                })
                .or_insert_with(|| {
                    // In this case, the hash was not present yet.
                    // we need to create the same color from small_colors
                    // into large_colors.
                    let ids = small_colors.indices(&color);
                    let new_color = large_colors.update(None, ids).unwrap();
                    assert_eq!(new_color, color);
                    new_color
                });
        });

        (large_hashes, large_colors)
    }

    pub fn search(
        &self,
        counter: SigCounter,
        similarity: bool,
        threshold: usize,
    ) -> Result<Vec<String>, Box<dyn std::error::Error>> {
        let mut matches = vec![];
        if similarity {
            unimplemented!("TODO: threshold correction")
        }

        for (dataset_id, size) in counter.most_common() {
            if size >= threshold {
                matches.push(self.sig_files[dataset_id as usize].to_str().unwrap().into());
            } else {
                break;
            };
        }
        Ok(matches)
    }

    pub fn gather(
        &self,
        mut counter: SigCounter,
        threshold: usize,
        query: &KmerMinHash,
    ) -> Result<Vec<GatherResult>, Box<dyn std::error::Error>> {
        let mut match_size = usize::max_value();
        let mut matches = vec![];

        while match_size > threshold && !counter.is_empty() {
            let (dataset_id, size) = counter.most_common()[0];
            match_size = if size >= threshold { size } else { break };

            let p;
            let match_path = if self.sig_files.is_empty() {
                p = PathBuf::new(); // TODO: Fix somehow?
                &p
            } else {
                &self.sig_files[dataset_id as usize]
            };

            let ref_match;
            let match_sig = if let Some(refsigs) = &self.ref_sigs {
                &refsigs[dataset_id as usize]
            } else {
                // TODO: remove swap_remove
                ref_match = Signature::from_path(&match_path)?.swap_remove(0);
                &ref_match
            };

            let mut match_mh = None;
            if let Some(Sketch::MinHash(mh)) = match_sig.select_sketch(&self.template) {
                match_mh = Some(mh);
            }
            let match_mh = match_mh.expect("Couldn't find a compatible MinHash");

            // Calculate stats
            let f_orig_query = match_size as f64 / query.size() as f64;
            let f_match = match_size as f64 / match_mh.size() as f64;
            let filename = match_path.to_str().unwrap().into();
            let name = match_sig.name();
            let unique_intersect_bp = match_mh.scaled() as usize * match_size;
            let gather_result_rank = matches.len();

            let (intersect_orig, _) = match_mh.intersection_size(query)?;
            let intersect_bp = (match_mh.scaled() as u64 * intersect_orig) as usize;

            let f_unique_to_query = intersect_orig as f64 / query.size() as f64;
            let match_ = match_sig.clone();

            // TODO: all of these
            let f_unique_weighted = 0.;
            let average_abund = 0;
            let median_abund = 0;
            let std_abund = 0;
            let md5 = "".into();
            let f_match_orig = 0.;
            let remaining_bp = 0;

            let result = GatherResult {
                intersect_bp,
                f_orig_query,
                f_match,
                f_unique_to_query,
                f_unique_weighted,
                average_abund,
                median_abund,
                std_abund,
                filename,
                name,
                md5,
                match_,
                f_match_orig,
                unique_intersect_bp,
                gather_result_rank,
                remaining_bp,
            };
            matches.push(result);

            // Prepare counter for finding the next match by decrementing
            // all hashes found in the current match in other datasets
            for hash in match_mh.iter_mins() {
                if let Some(color) = self.hash_to_color.get(hash) {
                    for dataset in self.colors.indices(color) {
                        counter.entry(*dataset).and_modify(|e| {
                            if *e > 0 {
                                *e -= 1
                            }
                        });
                    }
                }
            }
            counter.remove(&dataset_id);
        }
        Ok(matches)
    }

    pub fn counter_for_query(&self, query: &KmerMinHash) -> SigCounter {
        query
            .iter_mins()
            .filter_map(|hash| self.hash_to_color.get(hash))
            .flat_map(|color| self.colors.indices(color))
            .cloned()
            .collect()
    }

    pub fn counter(&self) -> SigCounter {
        self.hash_to_color
            .iter()
            .flat_map(|(_, color)| self.colors.indices(color))
            .cloned()
            .collect()
    }

    pub fn template(&self) -> Sketch {
        self.template.clone()
    }
    // TODO: mh should be a sketch, or even a sig...
    pub(crate) fn find_signatures(
        &self,
        mh: &KmerMinHash,
        threshold: f64,
        containment: bool,
        _ignore_scaled: bool,
    ) -> Result<Vec<(f64, Signature, String)>, Error> {
        /*
        let template_mh = None;
        if let Sketch::MinHash(mh) = self.template {
            template_mh = Some(mh);
        };
        // TODO: throw error
        let template_mh = template_mh.unwrap();

        let tmp_mh;
        let mh = if template_mh.scaled() > mh.scaled() {
            // TODO: proper error here
            tmp_mh = mh.downsample_scaled(self.scaled)?;
            &tmp_mh
        } else {
            mh
        };

                if self.scaled < mh.scaled() && !ignore_scaled {
                        return Err(LcaDBError::ScaledMismatchError {
                                db: self.scaled,
                                query: mh.scaled(),
                        }
                        .into());
                }
        */

        // TODO: proper threshold calculation
        let threshold: usize = (threshold * (mh.size() as f64)) as _;

        let counter = self.counter_for_query(&mh);

        debug!(
            "number of matching signatures for hashes: {}",
            counter.len()
        );

        let mut results = vec![];
        for (dataset_id, size) in counter.most_common() {
            let match_size = if size >= threshold { size } else { break };

            let p;
            let match_path = if self.sig_files.is_empty() {
                p = PathBuf::new(); // TODO: Fix somehow?
                &p
            } else {
                &self.sig_files[dataset_id as usize]
            };

            let ref_match;
            let match_sig = if let Some(refsigs) = &self.ref_sigs {
                &refsigs[dataset_id as usize]
            } else {
                // TODO: remove swap_remove
                ref_match = Signature::from_path(&match_path)?.swap_remove(0);
                &ref_match
            };

            let mut match_mh = None;
            if let Some(Sketch::MinHash(mh)) = match_sig.select_sketch(&self.template) {
                match_mh = Some(mh);
            }
            let match_mh = match_mh.unwrap();

            if size >= threshold {
                let score = if containment {
                    size as f64 / mh.size() as f64
                } else {
                    size as f64 / (mh.size() + match_size - size) as f64
                };
                let filename = match_path.to_str().unwrap().into();
                let mut sig = match_sig.clone();
                sig.reset_sketches();
                sig.push(Sketch::MinHash(match_mh.clone()));
                results.push((score, sig, filename));
            } else {
                break;
            };
        }
        Ok(results)
    }
}

#[derive(CopyGetters, Getters, Setters, Serialize, Deserialize, Debug)]
pub struct GatherResult {
    #[getset(get_copy = "pub")]
    intersect_bp: usize,

    #[getset(get_copy = "pub")]
    f_orig_query: f64,

    #[getset(get_copy = "pub")]
    f_match: f64,

    f_unique_to_query: f64,
    f_unique_weighted: f64,
    average_abund: usize,
    median_abund: usize,
    std_abund: usize,

    #[getset(get = "pub")]
    filename: String,

    #[getset(get = "pub")]
    name: String,

    md5: String,
    match_: Signature,
    f_match_orig: f64,
    unique_intersect_bp: usize,
    gather_result_rank: usize,
    remaining_bp: usize,
}

impl GatherResult {
    pub fn get_match(&self) -> Signature {
        self.match_.clone()
    }
}

impl<'a> Index<'a> for RevIndex {
    type Item = Signature;

    fn insert(&mut self, _node: Self::Item) -> Result<(), Error> {
        unimplemented!()
    }

    fn save<P: AsRef<Path>>(&self, _path: P) -> Result<(), Error> {
        unimplemented!()
    }

    fn load<P: AsRef<Path>>(_path: P) -> Result<(), Error> {
        unimplemented!()
    }

    fn signatures(&self) -> Vec<Self::Item> {
        if let Some(ref sigs) = self.ref_sigs {
            sigs.iter().map(|x| x.clone()).collect()
        } else {
            unimplemented!()
        }
    }

    fn signature_refs(&self) -> Vec<&Self::Item> {
        unimplemented!()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::sketch::minhash::max_hash_for_scaled;

    #[test]
    fn greyhound_new() {
        let max_hash = max_hash_for_scaled(10000);
        let template = Sketch::MinHash(
            KmerMinHash::builder()
                .num(0u32)
                .ksize(31)
                .max_hash(max_hash)
                .build(),
        );
        let search_sigs = [
            "../../tests/test-data/gather/GCF_000006945.2_ASM694v2_genomic.fna.gz.sig".into(),
            "../../tests/test-data/gather/GCF_000007545.1_ASM754v1_genomic.fna.gz.sig".into(),
        ];
        let index = RevIndex::new(&search_sigs, &template, 0, None, false);
        assert_eq!(index.colors.len(), 3);
    }

    #[test]
    fn greyhound_many() {
        let max_hash = max_hash_for_scaled(10000);
        let template = Sketch::MinHash(
            KmerMinHash::builder()
                .num(0u32)
                .ksize(31)
                .max_hash(max_hash)
                .build(),
        );
        let search_sigs = [
            "../../tests/test-data/gather/GCF_000006945.2_ASM694v2_genomic.fna.gz.sig".into(),
            "../../tests/test-data/gather/GCF_000007545.1_ASM754v1_genomic.fna.gz.sig".into(),
            "../../tests/test-data/gather/GCF_000008105.1_ASM810v1_genomic.fna.gz.sig".into(),
        ];

        let index = RevIndex::new(&search_sigs, &template, 0, None, false);
        /*
         dbg!(&index.colors.colors);
         0: 86
         1: 132
         2: 91
         (0, 1): 53
         (0, 2): 90
         (1, 2): 26
         (0, 1, 2): 261
         union: 739
        */
        //assert_eq!(index.colors.len(), 3);
        assert_eq!(index.colors.len(), 7);
    }
}
