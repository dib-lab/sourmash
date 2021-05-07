use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;
use std::fmt::Write;
use std::iter::{Iterator, Peekable};
use std::str;
use std::sync::Mutex;

use serde::de::Deserializer;
use serde::ser::{SerializeStruct, Serializer};
use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

use crate::_hash_murmur;
use crate::encodings::HashFunctions;
use crate::signature::SigsTrait;
use crate::sketch::hyperloglog::HyperLogLog;
use crate::Error;

#[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))]
use wasm_bindgen::prelude::*;

pub fn max_hash_for_scaled(scaled: u64) -> u64 {
    match scaled {
        0 => 0,
        1 => u64::max_value(),
        _ => (u64::max_value() as f64 / scaled as f64) as u64,
    }
}

pub fn scaled_for_max_hash(max_hash: u64) -> u64 {
    match max_hash {
        0 => 0,
        _ => (u64::max_value() as f64 / max_hash as f64) as u64,
    }
}

#[cfg_attr(all(target_arch = "wasm32", target_vendor = "unknown"), wasm_bindgen)]
#[derive(Debug, TypedBuilder)]
pub struct KmerMinHash {
    num: u32,
    ksize: u32,

    #[builder(setter(into), default = HashFunctions::murmur64_DNA)]
    hash_function: HashFunctions,

    #[builder(default = 42u64)]
    seed: u64,

    #[builder(default = u64::max_value())]
    max_hash: u64,

    #[builder(default)]
    mins: Vec<u64>,

    #[builder(default)]
    abunds: Option<Vec<u64>>,

    #[builder(default)]
    md5sum: Mutex<Option<String>>,
}

impl PartialEq for KmerMinHash {
    fn eq(&self, other: &KmerMinHash) -> bool {
        // TODO: check all other fields?
        self.md5sum() == other.md5sum()
    }
}

impl Clone for KmerMinHash {
    fn clone(&self) -> Self {
        KmerMinHash {
            num: self.num,
            ksize: self.ksize,
            hash_function: self.hash_function,
            seed: self.seed,
            max_hash: self.max_hash,
            mins: self.mins.clone(),
            abunds: self.abunds.clone(),
            md5sum: Mutex::new(Some(self.md5sum())),
        }
    }
}

impl Default for KmerMinHash {
    fn default() -> KmerMinHash {
        KmerMinHash {
            num: 1000,
            ksize: 21,
            hash_function: HashFunctions::murmur64_DNA,
            seed: 42,
            max_hash: 0,
            mins: Vec::with_capacity(1000),
            abunds: None,
            md5sum: Mutex::new(None),
        }
    }
}

impl Serialize for KmerMinHash {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let n_fields = match &self.abunds {
            Some(_) => 8,
            _ => 7,
        };

        let mut partial = serializer.serialize_struct("KmerMinHash", n_fields)?;
        partial.serialize_field("num", &self.num)?;
        partial.serialize_field("ksize", &self.ksize)?;
        partial.serialize_field("seed", &self.seed)?;
        partial.serialize_field("max_hash", &self.max_hash)?;
        partial.serialize_field("mins", &self.mins)?;
        partial.serialize_field("md5sum", &self.md5sum())?;

        if let Some(abunds) = &self.abunds {
            partial.serialize_field("abundances", abunds)?;
        }

        partial.serialize_field("molecule", &self.hash_function.to_string())?;

        partial.end()
    }
}

impl<'de> Deserialize<'de> for KmerMinHash {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct TempSig {
            num: u32,
            ksize: u32,
            seed: u64,
            max_hash: u64,
            md5sum: String,
            mins: Vec<u64>,
            abundances: Option<Vec<u64>>,
            molecule: String,
        }

        let tmpsig = TempSig::deserialize(deserializer)?;

        let num = if tmpsig.max_hash != 0 { 0 } else { tmpsig.num };
        let hash_function = match tmpsig.molecule.to_lowercase().as_ref() {
            "protein" => HashFunctions::murmur64_protein,
            "dayhoff" => HashFunctions::murmur64_dayhoff,
            "hp" => HashFunctions::murmur64_hp,
            "dna" => HashFunctions::murmur64_DNA,
            _ => unimplemented!(), // TODO: throw error here
        };

        // This shouldn't be necessary, but at some point we
        // created signatures with unordered mins =(
        let (mins, abunds) = if let Some(abunds) = tmpsig.abundances {
            let mut values: Vec<(_, _)> = tmpsig.mins.iter().zip(abunds.iter()).collect();
            values.sort();
            let mins = values.iter().map(|(v, _)| **v).collect();
            let abunds = values.iter().map(|(_, v)| **v).collect();
            (mins, Some(abunds))
        } else {
            let mut values: Vec<_> = tmpsig.mins.into_iter().collect();
            values.sort_unstable();
            (values, None)
        };

        Ok(KmerMinHash {
            num,
            ksize: tmpsig.ksize,
            seed: tmpsig.seed,
            max_hash: tmpsig.max_hash,
            md5sum: Mutex::new(Some(tmpsig.md5sum)),
            mins,
            abunds,
            hash_function,
        })
    }
}

impl KmerMinHash {
    pub fn new(
        scaled: u64,
        ksize: u32,
        hash_function: HashFunctions,
        seed: u64,
        track_abundance: bool,
        num: u32,
    ) -> KmerMinHash {
        let mins: Vec<u64>;
        let abunds: Option<Vec<u64>>;

        if num > 0 {
            mins = Vec::with_capacity(num as usize);
        } else {
            mins = Vec::with_capacity(1000);
        }

        if track_abundance {
            abunds = Some(Vec::with_capacity(mins.capacity()));
        } else {
            abunds = None
        }

        let max_hash = max_hash_for_scaled(scaled);

        KmerMinHash {
            num,
            ksize,
            hash_function,
            seed,
            max_hash,
            mins,
            abunds,
            md5sum: Mutex::new(None),
        }
    }

    pub fn num(&self) -> u32 {
        self.num
    }

    pub fn is_protein(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_protein
    }

    pub fn max_hash(&self) -> u64 {
        self.max_hash
    }

    pub fn scaled(&self) -> u64 {
        scaled_for_max_hash(self.max_hash)
    }

    pub fn clear(&mut self) {
        self.mins.clear();
        if let Some(ref mut abunds) = self.abunds {
            abunds.clear();
        }
    }

    pub fn is_empty(&self) -> bool {
        self.mins.is_empty()
    }

    pub fn set_hash_function(&mut self, h: HashFunctions) -> Result<(), Error> {
        if self.hash_function == h {
            return Ok(());
        }

        if !self.is_empty() {
            return Err(Error::NonEmptyMinHash {
                message: "hash_function".into(),
            });
        }

        self.hash_function = h;
        Ok(())
    }

    pub fn track_abundance(&self) -> bool {
        self.abunds.is_some()
    }

    pub fn enable_abundance(&mut self) -> Result<(), Error> {
        if !self.mins.is_empty() {
            return Err(Error::NonEmptyMinHash {
                message: "track_abundance=True".into(),
            });
        }

        self.abunds = Some(vec![]);

        Ok(())
    }

    pub fn disable_abundance(&mut self) {
        self.abunds = None;
    }

    fn reset_md5sum(&self) {
        let mut data = self.md5sum.lock().unwrap();
        if data.is_some() {
            *data = None;
        }
    }

    pub fn md5sum(&self) -> String {
        let mut data = self.md5sum.lock().unwrap();
        if data.is_none() {
            let mut buffer = String::with_capacity(20);

            let mut md5_ctx = md5::Context::new();
            write!(&mut buffer, "{}", self.ksize()).unwrap();
            md5_ctx.consume(&buffer);
            buffer.clear();
            for x in &self.mins {
                write!(&mut buffer, "{}", x).unwrap();
                md5_ctx.consume(&buffer);
                buffer.clear();
            }
            *data = Some(format!("{:x}", md5_ctx.compute()));
        }
        data.clone().unwrap()
    }

    pub fn add_hash(&mut self, hash: u64) {
        self.add_hash_with_abundance(hash, 1);
    }

    pub fn add_hash_with_abundance(&mut self, hash: u64, abundance: u64) {
        let current_max = match self.mins.last() {
            Some(&x) => x,
            None => u64::max_value(),
        };

        if hash > self.max_hash && self.max_hash != 0 {
            // This is a scaled minhash, and we don't need to add the new hash
            return;
        }

        if self.num == 0 && self.max_hash == 0 {
            // why did you create this minhash? it will always be empty...
            return;
        }

        if abundance == 0 {
            // well, don't add it.
            return;
        }

        // From this point on, hash is within scaled (or no scaled specified).

        // empty mins? add it.
        if self.mins.is_empty() {
            self.mins.push(hash);
            if let Some(ref mut abunds) = self.abunds {
                abunds.push(abundance);
                self.reset_md5sum();
            }
            return;
        }

        if hash <= self.max_hash || hash <= current_max || (self.mins.len() as u32) < self.num {
            // "good" hash - within range, smaller than current entry, or
            // still have space available
            let pos = match self.mins.binary_search(&hash) {
                Ok(p) => p,
                Err(p) => p,
            };

            if pos == self.mins.len() {
                // at end - must still be growing, we know the list won't
                // get too long
                self.mins.push(hash);
                self.reset_md5sum();
                if let Some(ref mut abunds) = self.abunds {
                    abunds.push(abundance);
                }
            } else if self.mins[pos] != hash {
                // didn't find hash in mins, so inserting somewhere
                // in the middle; shrink list if needed.
                self.mins.insert(pos, hash);
                if let Some(ref mut abunds) = self.abunds {
                    abunds.insert(pos, abundance);
                }

                // is it too big now?
                if self.num != 0 && self.mins.len() > (self.num as usize) {
                    self.mins.pop();
                    if let Some(ref mut abunds) = self.abunds {
                        abunds.pop();
                    }
                }
                self.reset_md5sum();
            } else if let Some(ref mut abunds) = self.abunds {
                // pos == hash: hash value already in mins, inc count by abundance
                abunds[pos] += abundance;
            }
        }
    }

    pub fn set_hash_with_abundance(&mut self, hash: u64, abundance: u64) {
        let mut found = false;
        if let Ok(pos) = self.mins.binary_search(&hash) {
            if self.mins[pos] == hash {
                found = true;
                if let Some(ref mut abunds) = self.abunds {
                    abunds[pos] = abundance;
                }
            }
        }

        if !found {
            self.add_hash_with_abundance(hash, abundance);
        }
    }

    pub fn add_word(&mut self, word: &[u8]) {
        let hash = _hash_murmur(word, self.seed);
        self.add_hash(hash);
    }

    pub fn remove_hash(&mut self, hash: u64) {
        if let Ok(pos) = self.mins.binary_search(&hash) {
            if self.mins[pos] == hash {
                self.mins.remove(pos);
                self.reset_md5sum();
                if let Some(ref mut abunds) = self.abunds {
                    abunds.remove(pos);
                }
            }
        };
    }

    pub fn remove_many(&mut self, hashes: &[u64]) -> Result<(), Error> {
        for min in hashes {
            self.remove_hash(*min);
        }
        Ok(())
    }

    pub fn merge(&mut self, other: &KmerMinHash) -> Result<(), Error> {
        self.check_compatible(other)?;
        let max_size = self.mins.len() + other.mins.len();
        let mut merged: Vec<u64> = Vec::with_capacity(max_size);
        let mut merged_abunds: Vec<u64> = Vec::with_capacity(max_size);

        {
            let mut self_iter = self.mins.iter();
            let mut other_iter = other.mins.iter();

            let mut self_abunds_iter: Option<std::slice::Iter<'_, u64>>;
            if let Some(ref mut abunds) = self.abunds {
                self_abunds_iter = Some(abunds.iter());
            } else {
                self_abunds_iter = None;
            }

            let mut other_abunds_iter: Option<std::slice::Iter<'_, u64>>;
            if let Some(ref abunds) = other.abunds {
                other_abunds_iter = Some(abunds.iter());
            } else {
                other_abunds_iter = None;
            }

            let mut self_value = self_iter.next();
            let mut other_value = other_iter.next();
            while self_value.is_some() {
                let value = self_value.unwrap();
                match other_value {
                    None => {
                        merged.push(*value);
                        merged.extend(self_iter);
                        if let Some(sai) = self_abunds_iter {
                            merged_abunds.extend(sai);
                        }
                        break;
                    }
                    Some(x) if x < value => {
                        merged.push(*x);
                        other_value = other_iter.next();

                        if let Some(ref mut oai) = other_abunds_iter {
                            if let Some(v) = oai.next() {
                                merged_abunds.push(*v)
                            }
                        }
                    }
                    Some(x) if x == value => {
                        merged.push(*x);
                        other_value = other_iter.next();
                        self_value = self_iter.next();

                        if let Some(ref mut oai) = other_abunds_iter {
                            if let Some(v) = oai.next() {
                                if let Some(ref mut sai) = self_abunds_iter {
                                    if let Some(s) = sai.next() {
                                        merged_abunds.push(*v + *s)
                                    }
                                }
                            }
                        }
                    }
                    Some(x) if x > value => {
                        merged.push(*value);
                        self_value = self_iter.next();

                        if let Some(ref mut sai) = self_abunds_iter {
                            if let Some(v) = sai.next() {
                                merged_abunds.push(*v)
                            }
                        }
                    }
                    Some(_) => {}
                }
            }
            if let Some(value) = other_value {
                merged.push(*value);
            }
            merged.extend(other_iter);
            if let Some(oai) = other_abunds_iter {
                merged_abunds.extend(oai);
            }
        }

        if merged.len() < (self.num as usize) || (self.num as usize) == 0 {
            self.mins = merged;
            self.abunds = if merged_abunds.is_empty() {
                if self.abunds.is_some() {
                    Some(vec![])
                } else {
                    None
                }
            } else {
                Some(merged_abunds)
            };
        } else {
            self.mins = merged.into_iter().take(self.num as usize).collect();
            self.abunds = if merged_abunds.is_empty() {
                if self.abunds.is_some() {
                    Some(vec![])
                } else {
                    None
                }
            } else {
                Some(merged_abunds.into_iter().take(self.num as usize).collect())
            }
        }

        self.reset_md5sum();
        Ok(())
    }

    pub fn add_from(&mut self, other: &KmerMinHash) -> Result<(), Error> {
        for min in &other.mins {
            self.add_hash(*min);
        }
        Ok(())
    }

    pub fn add_many(&mut self, hashes: &[u64]) -> Result<(), Error> {
        for min in hashes {
            self.add_hash(*min);
        }
        Ok(())
    }

    pub fn add_many_with_abund(&mut self, hashes: &[(u64, u64)]) -> Result<(), Error> {
        for item in hashes {
            self.add_hash_with_abundance(item.0, item.1);
        }
        Ok(())
    }

    pub fn count_common(&self, other: &KmerMinHash, downsample: bool) -> Result<u64, Error> {
        if downsample && self.max_hash != other.max_hash {
            let (first, second) = if self.max_hash < other.max_hash {
                (self, other)
            } else {
                (other, self)
            };
            let downsampled_mh = second.downsample_max_hash(first.max_hash)?;
            first.count_common(&downsampled_mh, false)
        } else {
            self.check_compatible(other)?;
            let iter = if self.size() < other.size() {
                Intersection::new(self.mins.iter(), other.mins.iter())
            } else {
                Intersection::new(other.mins.iter(), self.mins.iter())
            };

            Ok(iter.count() as u64)
        }
    }

    pub fn intersection(&self, other: &KmerMinHash) -> Result<(Vec<u64>, u64), Error> {
        self.check_compatible(other)?;

        let mut combined_mh = KmerMinHash::new(
            self.scaled(),
            self.ksize,
            self.hash_function,
            self.seed,
            self.abunds.is_some(),
            self.num,
        );

        combined_mh.merge(&self)?;
        combined_mh.merge(&other)?;

        let it1 = Intersection::new(self.mins.iter(), other.mins.iter());

        // TODO: there is probably a way to avoid this Vec here,
        // and pass the it1 as left in it2.
        let i1: Vec<u64> = it1.cloned().collect();
        let it2 = Intersection::new(i1.iter(), combined_mh.mins.iter());

        let common: Vec<u64> = it2.cloned().collect();
        Ok((common, combined_mh.mins.len() as u64))
    }

    // FIXME: intersection_size and count_common should be the same?
    // (for scaled minhashes)
    pub fn intersection_size(&self, other: &KmerMinHash) -> Result<(u64, u64), Error> {
        self.check_compatible(other)?;

        let mut combined_mh = KmerMinHash::new(
            self.scaled(),
            self.ksize,
            self.hash_function,
            self.seed,
            self.abunds.is_some(),
            self.num,
        );

        combined_mh.merge(&self)?;
        combined_mh.merge(&other)?;

        let it1 = Intersection::new(self.mins.iter(), other.mins.iter());

        // TODO: there is probably a way to avoid this Vec here,
        // and pass the it1 as left in it2.
        let i1: Vec<u64> = it1.cloned().collect();
        let it2 = Intersection::new(i1.iter(), combined_mh.mins.iter());

        Ok((it2.count() as u64, combined_mh.mins.len() as u64))
    }

    // calculate Jaccard similarity, ignoring abundance.
    pub fn jaccard(&self, other: &KmerMinHash) -> Result<f64, Error> {
        self.check_compatible(other)?;
        if let Ok((common, size)) = self.intersection_size(other) {
            Ok(common as f64 / u64::max(1, size) as f64)
        } else {
            Ok(0.0)
        }
    }

    // compare two minhashes, with abundance;
    // calculate their angular similarity.
    pub fn angular_similarity(&self, other: &KmerMinHash) -> Result<f64, Error> {
        self.check_compatible(other)?;

        if self.abunds.is_none() || other.abunds.is_none() {
            // TODO: throw error, we need abundance for this
            unimplemented!() // @CTB fixme
        }

        // TODO: check which one is smaller, swap around if needed

        let abunds = self.abunds.as_ref().unwrap();
        let other_abunds = other.abunds.as_ref().unwrap();

        let mut prod = 0;
        let mut other_iter = other.mins.iter().enumerate();
        let mut next_hash = other_iter.next();
        let a_sq: u64 = abunds.iter().map(|a| (a * a)).sum();
        let b_sq: u64 = other_abunds.iter().map(|a| (a * a)).sum();

        for (i, hash) in self.mins.iter().enumerate() {
            while let Some((j, k)) = next_hash {
                match k.cmp(hash) {
                    Ordering::Less => next_hash = other_iter.next(),
                    Ordering::Equal => {
                        // Calling `get_unchecked` here is safe since
                        // both `i` and `j` are valid indices
                        // (`i` and `j` came from valid iterator calls)
                        unsafe {
                            prod += abunds.get_unchecked(i) * other_abunds.get_unchecked(j);
                        }
                        break;
                    }
                    Ordering::Greater => break,
                }
            }
        }

        let norm_a = (a_sq as f64).sqrt();
        let norm_b = (b_sq as f64).sqrt();

        if norm_a == 0. || norm_b == 0. {
            return Ok(0.0);
        }
        let prod = f64::min(prod as f64 / (norm_a * norm_b), 1.);
        let distance = 2. * prod.acos() / PI;
        Ok(1. - distance)
    }

    pub fn similarity(
        &self,
        other: &KmerMinHash,
        ignore_abundance: bool,
        downsample: bool,
    ) -> Result<f64, Error> {
        if downsample && self.max_hash != other.max_hash {
            let (first, second) = if self.max_hash < other.max_hash {
                (self, other)
            } else {
                (other, self)
            };
            let downsampled_mh = second.downsample_max_hash(first.max_hash)?;
            first.similarity(&downsampled_mh, ignore_abundance, false)
        } else if ignore_abundance || self.abunds.is_none() || other.abunds.is_none() {
            self.jaccard(&other)
        } else {
            self.angular_similarity(&other)
        }
    }

    pub fn dayhoff(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_dayhoff
    }

    pub fn hp(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_hp
    }

    pub fn mins(&self) -> Vec<u64> {
        self.mins.clone()
    }

    pub fn iter_mins(&self) -> impl Iterator<Item = &u64> {
        self.mins.iter()
    }

    pub fn abunds(&self) -> Option<Vec<u64>> {
        self.abunds.clone()
    }

    // create a downsampled copy of self
    pub fn downsample_max_hash(&self, max_hash: u64) -> Result<KmerMinHash, Error> {
        let scaled = scaled_for_max_hash(max_hash);

        let mut new_mh = KmerMinHash::new(
            scaled,
            self.ksize,
            self.hash_function,
            self.seed,
            self.abunds.is_some(),
            self.num,
        );
        if self.abunds.is_some() {
            new_mh.add_many_with_abund(&self.to_vec_abunds())?;
        } else {
            new_mh.add_many(&self.mins)?;
        }
        Ok(new_mh)
    }

    pub fn to_vec_abunds(&self) -> Vec<(u64, u64)> {
        if let Some(abunds) = &self.abunds {
            self.mins
                .iter()
                .cloned()
                .zip(abunds.iter().cloned())
                .collect()
        } else {
            self.mins
                .iter()
                .cloned()
                .zip(std::iter::repeat(1))
                .collect()
        }
    }

    pub fn as_hll(&self) -> HyperLogLog {
        let mut hll = HyperLogLog::with_error_rate(0.01, self.ksize()).unwrap();

        for h in &self.mins {
            hll.add_hash(*h)
        }

        hll
    }

    // create a downsampled copy of self
    pub fn downsample_scaled(&self, scaled: u64) -> Result<KmerMinHash, Error> {
        let max_hash = max_hash_for_scaled(scaled);

        let mut new_mh = KmerMinHash::new(
            max_hash, // old max_hash => max_hash arg
            self.ksize,
            self.hash_function,
            self.seed,
            self.abunds.is_some(),
            self.num,
        );
        if self.abunds.is_some() {
            new_mh.add_many_with_abund(&self.to_vec_abunds())?;
        } else {
            new_mh.add_many(&self.mins)?;
        }
        Ok(new_mh)
    }
}

impl SigsTrait for KmerMinHash {
    fn size(&self) -> usize {
        self.mins.len()
    }

    fn to_vec(&self) -> Vec<u64> {
        self.mins.clone()
    }

    fn ksize(&self) -> usize {
        self.ksize as usize
    }

    fn seed(&self) -> u64 {
        self.seed
    }

    fn hash_function(&self) -> HashFunctions {
        self.hash_function
    }

    fn add_hash(&mut self, hash: u64) {
        self.add_hash_with_abundance(hash, 1);
    }

    fn check_compatible(&self, other: &KmerMinHash) -> Result<(), Error> {
        /*
        if self.num != other.num {
            return Err(Error::MismatchNum {
                n1: self.num,
                n2: other.num,
            }
            .into());
        }
        */
        if self.ksize != other.ksize {
            return Err(Error::MismatchKSizes);
        }
        if self.hash_function != other.hash_function {
            // TODO: fix this error
            return Err(Error::MismatchDNAProt);
        }
        if self.max_hash != other.max_hash {
            return Err(Error::MismatchScaled);
        }
        if self.seed != other.seed {
            return Err(Error::MismatchSeed);
        }
        Ok(())
    }
}

struct Intersection<T, I: Iterator<Item = T>> {
    iter: Peekable<I>,
    other: Peekable<I>,
}

impl<T, I: Iterator<Item = T>> Intersection<T, I> {
    pub fn new(left: I, right: I) -> Self {
        Intersection {
            iter: left.peekable(),
            other: right.peekable(),
        }
    }
}

impl<T: Ord, I: Iterator<Item = T>> Iterator for Intersection<T, I> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        loop {
            let res = match (self.iter.peek(), self.other.peek()) {
                (Some(ref left_key), Some(ref right_key)) => left_key.cmp(right_key),
                _ => return None,
            };

            match res {
                Ordering::Less => {
                    self.iter.next();
                }
                Ordering::Greater => {
                    self.other.next();
                }
                Ordering::Equal => {
                    self.other.next();
                    return self.iter.next();
                }
            }
        }
    }
}

//#############
// A MinHash implementation for low scaled or large cardinalities

#[cfg_attr(all(target_arch = "wasm32", target_vendor = "unknown"), wasm_bindgen)]
#[derive(Debug, TypedBuilder)]
pub struct KmerMinHashBTree {
    num: u32,
    ksize: u32,

    #[builder(setter(into), default = HashFunctions::murmur64_DNA)]
    hash_function: HashFunctions,

    #[builder(default = 42u64)]
    seed: u64,

    #[builder(default = u64::max_value())]
    max_hash: u64,

    #[builder(default)]
    mins: BTreeSet<u64>,

    #[builder(default)]
    abunds: Option<BTreeMap<u64, u64>>,

    #[builder(default = 0u64)]
    current_max: u64,

    #[builder(default)]
    md5sum: Mutex<Option<String>>,
}

impl PartialEq for KmerMinHashBTree {
    fn eq(&self, other: &KmerMinHashBTree) -> bool {
        // TODO: check all other fields?
        self.md5sum() == other.md5sum()
    }
}

impl Clone for KmerMinHashBTree {
    fn clone(&self) -> Self {
        KmerMinHashBTree {
            num: self.num,
            ksize: self.ksize,
            hash_function: self.hash_function,
            seed: self.seed,
            max_hash: self.max_hash,
            mins: self.mins.clone(),
            abunds: self.abunds.clone(),
            current_max: self.current_max,
            md5sum: Mutex::new(Some(self.md5sum())),
        }
    }
}

impl Default for KmerMinHashBTree {
    fn default() -> KmerMinHashBTree {
        KmerMinHashBTree {
            num: 1000,
            ksize: 21,
            hash_function: HashFunctions::murmur64_DNA,
            seed: 42,
            max_hash: 0,
            mins: Default::default(),
            abunds: None,
            current_max: 0,
            md5sum: Mutex::new(None),
        }
    }
}

impl Serialize for KmerMinHashBTree {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let n_fields = match &self.abunds {
            Some(_) => 8,
            _ => 7,
        };

        let mut partial = serializer.serialize_struct("KmerMinHashBTree", n_fields)?;
        partial.serialize_field("num", &self.num)?;
        partial.serialize_field("ksize", &self.ksize)?;
        partial.serialize_field("seed", &self.seed)?;
        partial.serialize_field("max_hash", &self.max_hash)?;
        partial.serialize_field("mins", &self.mins)?;
        partial.serialize_field("md5sum", &self.md5sum())?;

        if let Some(abunds) = &self.abunds {
            let abs: Vec<u64> = abunds.values().cloned().collect();
            partial.serialize_field("abundances", &abs)?;
        }

        partial.serialize_field("molecule", &self.hash_function.to_string())?;

        partial.end()
    }
}

impl<'de> Deserialize<'de> for KmerMinHashBTree {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct TempSig {
            num: u32,
            ksize: u32,
            seed: u64,
            max_hash: u64,
            md5sum: String,
            mins: Vec<u64>,
            abundances: Option<Vec<u64>>,
            molecule: String,
        }

        let tmpsig = TempSig::deserialize(deserializer)?;

        let num = if tmpsig.max_hash != 0 { 0 } else { tmpsig.num };
        let hash_function = match tmpsig.molecule.to_lowercase().as_ref() {
            "protein" => HashFunctions::murmur64_protein,
            "dayhoff" => HashFunctions::murmur64_dayhoff,
            "hp" => HashFunctions::murmur64_hp,
            "dna" => HashFunctions::murmur64_DNA,
            _ => unimplemented!(), // TODO: throw error here
        };

        let current_max;
        // This shouldn't be necessary, but at some point we
        // created signatures with unordered mins =(
        let (mins, abunds) = if let Some(abunds) = tmpsig.abundances {
            let mut values: Vec<(_, _)> = tmpsig.mins.iter().zip(abunds.iter()).collect();
            values.sort();
            let mins: BTreeSet<_> = values.iter().map(|(v, _)| **v).collect();
            let abunds = values.into_iter().map(|(v, x)| (*v, *x)).collect();
            current_max = *mins.iter().rev().next().unwrap_or(&0);
            (mins, Some(abunds))
        } else {
            current_max = 0;
            (tmpsig.mins.into_iter().collect(), None)
        };

        Ok(KmerMinHashBTree {
            num,
            ksize: tmpsig.ksize,
            seed: tmpsig.seed,
            max_hash: tmpsig.max_hash,
            md5sum: Mutex::new(Some(tmpsig.md5sum)),
            mins,
            abunds,
            hash_function,
            current_max,
        })
    }
}

impl KmerMinHashBTree {
    pub fn new(
        scaled: u64,
        ksize: u32,
        hash_function: HashFunctions,
        seed: u64,
        track_abundance: bool,
        num: u32,
    ) -> KmerMinHashBTree {
        let mins = Default::default();

        let abunds = if track_abundance {
            Some(Default::default())
        } else {
            None
        };

        let max_hash = max_hash_for_scaled(scaled);

        KmerMinHashBTree {
            num,
            ksize,
            hash_function,
            seed,
            max_hash,
            mins,
            abunds,
            current_max: 0,
            md5sum: Mutex::new(None),
        }
    }

    pub fn num(&self) -> u32 {
        self.num
    }

    pub fn is_protein(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_protein
    }

    pub fn max_hash(&self) -> u64 {
        self.max_hash
    }

    pub fn scaled(&self) -> u64 {
        scaled_for_max_hash(self.max_hash)
    }

    pub fn clear(&mut self) {
        self.mins.clear();
        if let Some(ref mut abunds) = self.abunds {
            abunds.clear();
        }
        self.current_max = 0;
    }

    pub fn is_empty(&self) -> bool {
        self.mins.is_empty()
    }

    pub fn set_hash_function(&mut self, h: HashFunctions) -> Result<(), Error> {
        if self.hash_function == h {
            return Ok(());
        }

        if !self.is_empty() {
            return Err(Error::NonEmptyMinHash {
                message: "hash_function".into(),
            });
        }

        self.hash_function = h;
        Ok(())
    }

    pub fn track_abundance(&self) -> bool {
        self.abunds.is_some()
    }

    pub fn enable_abundance(&mut self) -> Result<(), Error> {
        if !self.mins.is_empty() {
            return Err(Error::NonEmptyMinHash {
                message: "track_abundance=True".into(),
            });
        }

        self.abunds = Some(Default::default());

        Ok(())
    }

    pub fn disable_abundance(&mut self) {
        self.abunds = None;
    }

    fn reset_md5sum(&self) {
        let mut data = self.md5sum.lock().unwrap();
        if data.is_some() {
            *data = None;
        }
    }

    pub fn md5sum(&self) -> String {
        let mut data = self.md5sum.lock().unwrap();
        if data.is_none() {
            let mut buffer = String::with_capacity(20);

            let mut md5_ctx = md5::Context::new();
            write!(&mut buffer, "{}", self.ksize()).unwrap();
            md5_ctx.consume(&buffer);
            buffer.clear();
            for x in &self.mins {
                write!(&mut buffer, "{}", x).unwrap();
                md5_ctx.consume(&buffer);
                buffer.clear();
            }
            *data = Some(format!("{:x}", md5_ctx.compute()));
        }
        data.clone().unwrap()
    }

    pub fn add_hash_with_abundance(&mut self, hash: u64, abundance: u64) {
        if hash > self.max_hash && self.max_hash != 0 {
            // This is a scaled minhash, and we don't need to add the new hash
            return;
        }

        if self.num == 0 && self.max_hash == 0 {
            // why did you create this minhash? it will always be empty...
            return;
        }

        if abundance == 0 {
            // well, don't add it.
            return;
        }

        // From this point on, hash is within scaled (or no scaled specified).

        // empty mins? add it.
        if self.mins.is_empty() {
            self.mins.insert(hash);
            self.reset_md5sum();
            if let Some(ref mut abunds) = self.abunds {
                abunds.insert(hash, abundance);
            }
            self.current_max = hash;
            return;
        }

        if hash <= self.max_hash || hash <= self.current_max || (self.mins.len() as u32) < self.num
        {
            // "good" hash - within range, smaller than current entry, or
            // still have space available
            if self.mins.insert(hash) {
                self.reset_md5sum();
                if hash > self.current_max {
                    self.current_max = hash;
                }
            }
            if let Some(ref mut abunds) = self.abunds {
                *abunds.entry(hash).or_insert(0) += abundance;
            }

            // is it too big now?
            if self.num != 0 && self.mins.len() > (self.num as usize) {
                let last = *self.mins.iter().rev().next().unwrap();
                self.mins.remove(&last);
                self.reset_md5sum();
                if let Some(ref mut abunds) = self.abunds {
                    abunds.remove(&last);
                }
                self.current_max = *self.mins.iter().rev().next().unwrap();
            }
        }
    }

    pub fn add_word(&mut self, word: &[u8]) {
        let hash = _hash_murmur(word, self.seed);
        self.add_hash(hash);
    }

    pub fn remove_hash(&mut self, hash: u64) {
        if self.mins.remove(&hash) {
            self.reset_md5sum();
            if let Some(ref mut abunds) = self.abunds {
                abunds.remove(&hash);
            }
        }
        if hash == self.current_max {
            self.current_max = *self.mins.iter().rev().next().unwrap_or(&0);
        }
    }

    pub fn remove_many(&mut self, hashes: &[u64]) -> Result<(), Error> {
        for min in hashes {
            self.remove_hash(*min);
        }
        Ok(())
    }

    pub fn merge(&mut self, other: &KmerMinHashBTree) -> Result<(), Error> {
        self.check_compatible(other)?;
        let union = self.mins.union(&other.mins);

        let to_take = if self.num == 0 {
            usize::max_value()
        } else {
            self.num as usize
        };

        self.mins = union.take(to_take).cloned().collect();

        if let Some(abunds) = &self.abunds {
            if let Some(oabunds) = &other.abunds {
                let mut new_abunds = BTreeMap::new();

                for hash in &self.mins {
                    *new_abunds.entry(*hash).or_insert(0) +=
                        abunds.get(&hash).unwrap_or(&0) + oabunds.get(&hash).unwrap_or(&0);
                }
                self.abunds = Some(new_abunds)
            }
        }
        // Better safe than sorry, but could check in other places to avoid
        // always resetting
        self.reset_md5sum();

        Ok(())
    }

    pub fn add_from(&mut self, other: &KmerMinHashBTree) -> Result<(), Error> {
        for min in &other.mins {
            self.add_hash(*min);
        }
        Ok(())
    }

    pub fn add_many(&mut self, hashes: &[u64]) -> Result<(), Error> {
        for min in hashes {
            self.add_hash(*min);
        }
        Ok(())
    }

    pub fn add_many_with_abund(&mut self, hashes: &[(u64, u64)]) -> Result<(), Error> {
        for item in hashes {
            self.add_hash_with_abundance(item.0, item.1);
        }
        Ok(())
    }

    pub fn count_common(&self, other: &KmerMinHashBTree, downsample: bool) -> Result<u64, Error> {
        if downsample && self.max_hash != other.max_hash {
            let (first, second) = if self.max_hash < other.max_hash {
                (self, other)
            } else {
                (other, self)
            };
            let downsampled_mh = second.downsample_max_hash(first.max_hash)?;
            first.count_common(&downsampled_mh, false)
        } else {
            self.check_compatible(other)?;
            let iter = if self.size() < other.size() {
                Intersection::new(self.mins.iter(), other.mins.iter())
            } else {
                Intersection::new(other.mins.iter(), self.mins.iter())
            };

            Ok(iter.count() as u64)
        }
    }

    pub fn intersection(&self, other: &KmerMinHashBTree) -> Result<(Vec<u64>, u64), Error> {
        self.check_compatible(other)?;

        let mut combined_mh = KmerMinHashBTree::new(
            self.scaled(),
            self.ksize,
            self.hash_function,
            self.seed,
            self.abunds.is_some(),
            self.num,
        );

        combined_mh.merge(&self)?;
        combined_mh.merge(&other)?;

        let it1 = Intersection::new(self.mins.iter(), other.mins.iter());

        // TODO: there is probably a way to avoid this Vec here,
        // and pass the it1 as left in it2.
        let i1: Vec<u64> = it1.cloned().collect();
        let i2: Vec<u64> = combined_mh.mins.iter().cloned().collect();
        let it2 = Intersection::new(i1.iter(), i2.iter());

        let common: Vec<u64> = it2.cloned().collect();
        Ok((common, combined_mh.mins.len() as u64))
    }

    pub fn intersection_size(&self, other: &KmerMinHashBTree) -> Result<(u64, u64), Error> {
        self.check_compatible(other)?;

        let mut combined_mh = KmerMinHashBTree::new(
            self.scaled(),
            self.ksize,
            self.hash_function,
            self.seed,
            self.abunds.is_some(),
            self.num,
        );

        combined_mh.merge(&self)?;
        combined_mh.merge(&other)?;

        let it1 = Intersection::new(self.mins.iter(), other.mins.iter());

        // TODO: there is probably a way to avoid this Vec here,
        // and pass the it1 as left in it2.
        let i1: Vec<u64> = it1.cloned().collect();
        let i2: Vec<u64> = combined_mh.mins.iter().cloned().collect();
        let it2 = Intersection::new(i1.iter(), i2.iter());

        Ok((it2.count() as u64, combined_mh.mins.len() as u64))
    }

    // calculate Jaccard similarity, ignoring abundance.
    pub fn jaccard(&self, other: &KmerMinHashBTree) -> Result<f64, Error> {
        self.check_compatible(other)?;
        if let Ok((common, size)) = self.intersection_size(other) {
            Ok(common as f64 / u64::max(1, size) as f64)
        } else {
            Ok(0.0)
        }
    }

    // compare two minhashes, with abundance;
    // calculate their angular similarity.
    pub fn angular_similarity(&self, other: &KmerMinHashBTree) -> Result<f64, Error> {
        self.check_compatible(other)?;

        if self.abunds.is_none() || other.abunds.is_none() {
            // TODO: throw error, we need abundance for this
            unimplemented!() // @CTB fixme
        }

        let abunds = self.abunds.as_ref().unwrap();
        let other_abunds = other.abunds.as_ref().unwrap();

        let mut prod = 0;
        let a_sq: u64 = abunds.values().map(|a| (a * a)).sum();
        let b_sq: u64 = other_abunds.values().map(|a| (a * a)).sum();

        for (hash, value) in abunds.iter() {
            if let Some(oa) = other_abunds.get(&hash) {
                prod += value * oa
            }
        }

        let norm_a = (a_sq as f64).sqrt();
        let norm_b = (b_sq as f64).sqrt();

        if norm_a == 0. || norm_b == 0. {
            return Ok(0.0);
        }
        let prod = f64::min(prod as f64 / (norm_a * norm_b), 1.);
        let distance = 2. * prod.acos() / PI;
        Ok(1. - distance)
    }

    pub fn similarity(
        &self,
        other: &KmerMinHashBTree,
        ignore_abundance: bool,
        downsample: bool,
    ) -> Result<f64, Error> {
        if downsample && self.max_hash != other.max_hash {
            let (first, second) = if self.max_hash < other.max_hash {
                (self, other)
            } else {
                (other, self)
            };
            let downsampled_mh = second.downsample_max_hash(first.max_hash)?;
            first.similarity(&downsampled_mh, ignore_abundance, false)
        } else if ignore_abundance || self.abunds.is_none() || other.abunds.is_none() {
            self.jaccard(&other)
        } else {
            self.angular_similarity(&other)
        }
    }

    pub fn dayhoff(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_dayhoff
    }

    pub fn hp(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_hp
    }

    pub fn hash_function(&self) -> HashFunctions {
        self.hash_function
    }

    pub fn mins(&self) -> Vec<u64> {
        self.mins.iter().cloned().collect()
    }

    pub fn iter_mins(&self) -> impl Iterator<Item = &u64> {
        self.mins.iter()
    }

    pub fn abunds(&self) -> Option<Vec<u64>> {
        if let Some(abunds) = &self.abunds {
            Some(abunds.values().cloned().collect())
        } else {
            None
        }
    }

    // create a downsampled copy of self
    pub fn downsample_max_hash(&self, max_hash: u64) -> Result<KmerMinHashBTree, Error> {
        let scaled = scaled_for_max_hash(max_hash);

        let mut new_mh = KmerMinHashBTree::new(
            scaled,
            self.ksize,
            self.hash_function,
            self.seed,
            self.abunds.is_some(),
            self.num,
        );
        if self.abunds.is_some() {
            new_mh.add_many_with_abund(&self.to_vec_abunds())?;
        } else {
            new_mh.add_many(&self.mins())?;
        }
        Ok(new_mh)
    }

    pub fn to_vec_abunds(&self) -> Vec<(u64, u64)> {
        if let Some(abunds) = &self.abunds {
            abunds.iter().map(|(a, b)| (*a, *b)).collect()
        } else {
            self.mins
                .iter()
                .cloned()
                .zip(std::iter::repeat(1))
                .collect()
        }
    }
}

impl SigsTrait for KmerMinHashBTree {
    fn size(&self) -> usize {
        self.mins.len()
    }

    fn to_vec(&self) -> Vec<u64> {
        self.mins()
    }

    fn ksize(&self) -> usize {
        self.ksize as usize
    }

    fn seed(&self) -> u64 {
        self.seed
    }

    fn hash_function(&self) -> HashFunctions {
        self.hash_function
    }

    fn add_hash(&mut self, hash: u64) {
        self.add_hash_with_abundance(hash, 1);
    }

    fn check_compatible(&self, other: &KmerMinHashBTree) -> Result<(), Error> {
        /*
        if self.num != other.num {
            return Err(Error::MismatchNum {
                n1: self.num,
                n2: other.num,
            }
            .into());
        }
        */
        if self.ksize != other.ksize {
            return Err(Error::MismatchKSizes);
        }
        if self.hash_function != other.hash_function {
            // TODO: fix this error
            return Err(Error::MismatchDNAProt);
        }
        if self.max_hash != other.max_hash {
            return Err(Error::MismatchScaled);
        }
        if self.seed != other.seed {
            return Err(Error::MismatchSeed);
        }
        Ok(())
    }
}

impl From<KmerMinHashBTree> for KmerMinHash {
    fn from(other: KmerMinHashBTree) -> KmerMinHash {
        let mut new_mh = KmerMinHash::new(
            other.scaled(),
            other.ksize() as u32,
            other.hash_function(),
            other.seed(),
            other.track_abundance(),
            other.num(),
        );

        let mins = other.mins.into_iter().collect();
        let abunds = if let Some(abunds) = other.abunds {
            Some(abunds.values().cloned().collect())
        } else {
            None
        };

        new_mh.mins = mins;
        new_mh.abunds = abunds;

        new_mh
    }
}

impl From<KmerMinHash> for KmerMinHashBTree {
    fn from(other: KmerMinHash) -> KmerMinHashBTree {
        let mut new_mh = KmerMinHashBTree::new(
            other.scaled(),
            other.ksize() as u32,
            other.hash_function(),
            other.seed(),
            other.track_abundance(),
            other.num(),
        );

        let mins: BTreeSet<u64> = other.mins.into_iter().collect();
        let abunds = if let Some(abunds) = other.abunds {
            Some(mins.iter().cloned().zip(abunds.into_iter()).collect())
        } else {
            None
        };

        new_mh.mins = mins;
        new_mh.abunds = abunds;

        new_mh
    }
}
