use std::f64::consts::PI;
use std::fs::File;
use std::io;
use std::path::Path;
use std::slice;

use byteorder::{BigEndian, ByteOrder, LittleEndian, ReadBytesExt, WriteBytesExt};
use fixedbitset::FixedBitSet;

use crate::index::sbt::Update;
use crate::index::Comparable;
use crate::sketch::minhash::KmerMinHash;
use crate::Error;
use crate::HashIntoType;

#[derive(Debug, Default, Clone)]
pub struct Nodegraph {
    bs: Vec<FixedBitSet>,
    ksize: usize,
    occupied_bins: usize,
    unique_kmers: usize,
}

// TODO: not checking for unique_kmers,
// since it is not saved in a khmer nodegraph
impl PartialEq for Nodegraph {
    fn eq(&self, other: &Nodegraph) -> bool {
        self.bs == other.bs
            && self.occupied_bins == other.occupied_bins
            && self.ksize == other.ksize
    }
}

impl Update<Nodegraph> for Nodegraph {
    fn update(&self, other: &mut Nodegraph) -> Result<(), Error> {
        other.occupied_bins = other
            .bs
            .iter_mut()
            .zip(&self.bs)
            .enumerate()
            .map(|(i, (bs, bs_me))| {
                bs.union_with(bs_me);
                if i == 0 {
                    bs.count_ones(..)
                } else {
                    0
                }
            })
            .sum();
        Ok(())
    }
}

impl Update<Nodegraph> for KmerMinHash {
    fn update(&self, other: &mut Nodegraph) -> Result<(), Error> {
        for h in self.mins() {
            other.count(h);
        }
        Ok(())
    }
}

impl Nodegraph {
    pub fn new(tablesizes: &[usize], ksize: usize) -> Nodegraph {
        let mut bs = Vec::with_capacity(tablesizes.len());
        for size in tablesizes.iter() {
            bs.push(FixedBitSet::with_capacity(*size));
        }

        Nodegraph {
            bs,
            ksize,
            occupied_bins: 0,
            unique_kmers: 0,
        }
    }

    pub fn with_tables(tablesize: usize, n_tables: usize, ksize: usize) -> Nodegraph {
        let mut tablesizes = Vec::with_capacity(n_tables);

        let mut i = u64::max((tablesize - 1) as u64, 2);
        if i % 2 == 0 {
            i -= 1
        }

        while tablesizes.len() != n_tables {
            if primal_check::miller_rabin(i) {
                tablesizes.push(i as usize);
            }
            if i == 1 {
                break;
            }
            i -= 2;
        }

        Nodegraph::new(tablesizes.as_slice(), ksize)
    }

    pub(crate) fn count_kmer(&mut self, kmer: &[u8]) -> bool {
        let h = _hash(kmer);
        self.count(h)
    }

    pub fn count(&mut self, hash: HashIntoType) -> bool {
        let mut is_new_kmer = false;

        for (i, bitset) in self.bs.iter_mut().enumerate() {
            let bin = hash % bitset.len() as u64;
            if !bitset.put(bin as usize) {
                if i == 0 {
                    self.occupied_bins += 1;
                }
                is_new_kmer = true;
            }
        }

        if is_new_kmer {
            self.unique_kmers += 1
        }
        is_new_kmer
    }

    pub fn get(&self, hash: HashIntoType) -> usize {
        for bitset in &self.bs {
            let bin = hash % bitset.len() as u64;
            if !bitset.contains(bin as usize) {
                return 0;
            }
        }
        1
    }

    pub(crate) fn get_kmer(&self, kmer: &[u8]) -> usize {
        let h = _hash(kmer);
        self.get(h)
    }

    pub fn expected_collisions(&self) -> f64 {
        let min_size = self.bs.iter().map(|x| x.len()).min().unwrap();
        let n_ht = self.bs.len();
        let occupancy = self.occupied_bins;

        let fp_one = occupancy as f64 / min_size as f64;
        f64::powf(fp_one as f64, n_ht as f64)
    }

    pub fn tablesize(&self) -> usize {
        self.bs.iter().map(|x| x.len()).sum()
    }

    pub fn noccupied(&self) -> usize {
        self.occupied_bins
    }

    pub fn matches(&self, mh: &KmerMinHash) -> usize {
        mh.iter_mins().filter(|x| self.get(**x) == 1).count()
    }

    /// upper-bound estimate of the angular similarity with an abundance MinHash
    ///
    /// Note that this is not a *tight* bound: it might overestimate a lot.
    /// The current goal is guaranteeing it doesn't *underestimate*
    /// (and, if possible, it doesn't overestimate all the time).
    pub fn angular_similarity_upper_bound(&self, other: &KmerMinHash) -> Result<f64, Error> {
        if !other.track_abundance() {
            // TODO: throw error, we need abundance for this
            unimplemented!()
        }
        let other_abunds = other.to_vec_abunds();

        let mut prod = 0;
        let mut a_sq = 0_u64;
        let b_sq: u64 = other_abunds.iter().map(|(_, b)| (b * b)).sum();

        for (hash, abund) in other_abunds {
            if self.get(hash) == 1 {
                // TODO: which one overestimate less?
                a_sq += abund * abund;
                //a_sq += 1;

                prod += abund * abund;
            }
        }

        let norm_a = (a_sq as f64).sqrt();
        let norm_b = (b_sq as f64).sqrt();

        if norm_a == 0. || norm_b == 0. {
            return Ok(0.0);
        }
        let prod = f64::min(prod as f64 / (norm_a * norm_b), 1.);
        let distance = 2. * prod.acos() / PI;
        // Adding some leeway with a small epsilon
        Ok(1. - distance + 1e-3)
    }

    pub fn ntables(&self) -> usize {
        self.bs.len()
    }

    pub fn ksize(&self) -> usize {
        self.ksize
    }

    pub fn into_bitsets(self) -> Vec<FixedBitSet> {
        self.bs
    }

    // save
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
        // TODO: if it ends with gz, open a compressed file
        // might use get_output here?
        self.save_to_writer(&mut File::create(path)?)?;
        Ok(())
    }

    pub fn save_to_writer<W>(&self, wtr: &mut W) -> Result<(), Error>
    where
        W: io::Write,
    {
        wtr.write_all(b"OXLI")?;
        wtr.write_u8(4)?; // version
        wtr.write_u8(2)?; // ht_type
        wtr.write_u32::<LittleEndian>(self.ksize as u32)?; // ksize
        wtr.write_u8(self.bs.len() as u8)?; // n_tables
        wtr.write_u64::<LittleEndian>(self.occupied_bins as u64)?; // n_occupied
        for count in &self.bs {
            let tablesize = count.len();
            wtr.write_u64::<LittleEndian>(tablesize as u64)?;

            let byte_size = tablesize / 8 + 1;
            let (div, rem) = (byte_size / 4, byte_size % 4);

            // Once this issue and PR are solved, this is a one liner:
            // https://github.com/BurntSushi/byteorder/issues/155
            // https://github.com/BurntSushi/byteorder/pull/166
            //wtr.write_u32_from::<LittleEndian>(&count.as_slice()[..div])?;
            let slice = &count.as_slice()[..div];
            let buf = unsafe {
                use std::mem::size_of;

                let len = size_of::<u32>() * slice.len();
                slice::from_raw_parts(slice.as_ptr() as *const u8, len)
            };
            wtr.write_all(&buf)?;
            // Replace when byteorder PR is released

            if rem != 0 {
                let mut cursor = [0u8; 4];
                LittleEndian::write_u32(&mut cursor, count.as_slice()[div]);
                for item in cursor.iter().take(rem) {
                    wtr.write_u8(*item)?;
                }
            }
        }
        Ok(())
    }

    pub fn from_reader<R>(rdr: R) -> Result<Nodegraph, Error>
    where
        R: io::Read,
    {
        let (mut rdr, _format) = niffler::get_reader(Box::new(rdr))?;

        let signature = rdr.read_u32::<BigEndian>()?;
        assert_eq!(signature, 0x4f58_4c49);

        let version = rdr.read_u8()?;
        assert_eq!(version, 0x04);

        let ht_type = rdr.read_u8()?;
        assert_eq!(ht_type, 0x02);

        let ksize = rdr.read_u32::<LittleEndian>()?;
        let n_tables = rdr.read_u8()?;
        let occupied_bins = rdr.read_u64::<LittleEndian>()? as usize;

        let mut bs = Vec::with_capacity(n_tables as usize);
        for _i in 0..n_tables {
            let tablesize: usize = rdr.read_u64::<LittleEndian>()? as usize;
            let byte_size = tablesize / 8 + 1;

            let rem = byte_size % 4;
            let blocks: Vec<u32> = if rem == 0 {
                let mut blocks = vec![0; byte_size / 4];
                rdr.read_u32_into::<LittleEndian>(&mut blocks)?;
                blocks
            } else {
                let mut blocks = vec![0; byte_size / 4];
                rdr.read_u32_into::<LittleEndian>(&mut blocks)?;

                let mut values = [0u8; 4];
                for item in values.iter_mut().take(rem) {
                    let byte = rdr.read_u8().expect("error reading bins");
                    *item = byte;
                }
                let mut block = vec![0u32; 1];
                LittleEndian::read_u32_into(&values, &mut block);
                blocks.push(block[0]);
                blocks
            };

            let counts = FixedBitSet::with_capacity_and_blocks(tablesize, blocks);
            bs.push(counts);
        }

        Ok(Nodegraph {
            bs,
            ksize: ksize as usize,
            occupied_bins,
            unique_kmers: 0, // This is a khmer issue, it doesn't save unique_kmers
        })
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Nodegraph, Error> {
        let mut reader = io::BufReader::new(File::open(path)?);
        Nodegraph::from_reader(&mut reader)
    }

    pub fn tablesizes(&self) -> Vec<u64> {
        self.bs.iter().map(|x| x.len() as u64).collect()
    }

    pub fn n_occupied_bins(&self) -> usize {
        self.occupied_bins
    }

    pub fn unique_kmers(&self) -> usize {
        self.unique_kmers
    }
}

impl Comparable<&Nodegraph> for Nodegraph {
    fn similarity(&self, other: &&Nodegraph) -> f64 {
        self.similarity(*other)
    }

    fn containment(&self, other: &&Nodegraph) -> f64 {
        self.containment(*other)
    }
}

impl Comparable<Nodegraph> for Nodegraph {
    fn similarity(&self, other: &Nodegraph) -> f64 {
        let result: usize = self
            .bs
            .iter()
            .zip(&other.bs)
            .map(|(bs, bs_other)| bs.intersection(bs_other).count())
            .sum();
        let size: usize = self
            .bs
            .iter()
            .zip(&other.bs)
            .map(|(bs, bs_other)| bs.union(bs_other).count())
            .sum();
        result as f64 / size as f64
    }

    fn containment(&self, other: &Nodegraph) -> f64 {
        let result: usize = self
            .bs
            .iter()
            .zip(&other.bs)
            .map(|(bs, bs_other)| bs.intersection(bs_other).count())
            .sum();
        let size: usize = self.bs.iter().map(|bs| bs.count_ones(..)).sum();
        result as f64 / size as f64
    }
}

impl Comparable<KmerMinHash> for Nodegraph {
    fn similarity(&self, _other: &KmerMinHash) -> f64 {
        unimplemented!()
    }

    fn containment(&self, _other: &KmerMinHash) -> f64 {
        /*
          let result: usize = other.mins().iter().map(|h| self.get(*h)).sum();
          result as f64 / self.size() as f64
        */
        unimplemented!()
    }
}

fn twobit_repr(a: u8) -> HashIntoType {
    match a as char {
        'A' => 0,
        'C' => 2,
        'G' => 3,
        'T' => 1,
        _ => unimplemented!(),
    }
}

fn twobit_comp(a: u8) -> HashIntoType {
    match a as char {
        'A' => 1,
        'C' => 3,
        'G' => 2,
        'T' => 0,
        _ => unimplemented!(),
    }
}

fn uniqify_rc(f: HashIntoType, r: HashIntoType) -> HashIntoType {
    if f < r {
        f
    } else {
        r
    }
}

fn _hash(kmer: &[u8]) -> HashIntoType {
    let ksize = kmer.len();
    let mut hash = 0;
    let mut rev = 0;

    hash |= twobit_repr(kmer[0]);
    rev |= twobit_comp(kmer[ksize - 1]);

    let mut i = 1;
    let mut j: isize = (ksize - 2) as isize;

    while i < ksize {
        hash <<= 2;
        rev <<= 2;

        hash |= twobit_repr(kmer[i]);
        rev |= twobit_comp(kmer[j as usize]);

        i += 1;
        j -= 1;
    }

    uniqify_rc(hash, rev)
}

#[cfg(test)]
mod test {
    use super::*;
    use std::io::{BufReader, BufWriter};
    use std::path::PathBuf;

    use proptest::collection::vec;
    use proptest::num::u64;
    use proptest::proptest;

    // Generate with khmer:
    // >>> a = khmer.Nodegraph(3, 23, 6)
    // >>> a.count("ACG")
    // >>> a.count("TTA")
    // >>> a.count("CGA")
    // >>> a.save("test.ng")
    // and dumping test.ng with xxd:
    // $ xxd -i test.ng
    static RAW_DATA: &[u8] = &[
        0x4f, 0x58, 0x4c, 0x49, 0x04, 0x02, 0x03, 0x00, 0x00, 0x00, 0x06, 0x03, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x13, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x09, 0x01,
        0x11, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0c, 0x01, 0x0d, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x0a, 0x08, 0x0b, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x21,
        0x00, 0x07, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x54, 0x05, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x06,
    ];

    static COMPRESSED_RAW_DATA: &[u8] = &[
        0x1f, 0x8b, 0x08, 0x08, 0x73, 0x88, 0x9f, 0x5e, 0x00, 0x03, 0x74, 0x65, 0x73, 0x74, 0x2e,
        0x6e, 0x67, 0x00, 0xf3, 0x8f, 0xf0, 0xf1, 0x64, 0x61, 0x62, 0x66, 0x60, 0x60, 0x60, 0x03,
        0x11, 0x20, 0x20, 0x0c, 0xa5, 0x19, 0x38, 0x19, 0x05, 0x61, 0x4c, 0x1e, 0x46, 0x5e, 0x28,
        0x8b, 0x8b, 0x83, 0x1b, 0xca, 0x52, 0x64, 0x60, 0x87, 0xb2, 0x42, 0x58, 0xa1, 0x0c, 0x36,
        0x00, 0x8d, 0xf0, 0xa9, 0x8b, 0x4f, 0x00, 0x00, 0x00,
    ];

    proptest! {
      #[test]
      fn count_and_get(hashes in vec(u64::ANY, 1..500)) {
          let mut ng: Nodegraph = Nodegraph::new(&[1000], 3);
          for hash in hashes {
              ng.count(hash);
              assert_eq!(ng.get(hash), 1);
          }
      }
    }

    #[test]
    fn load_compressed() {
        let mut reader = BufReader::new(&COMPRESSED_RAW_DATA[..]);

        let ng: Nodegraph = Nodegraph::from_reader(&mut reader).expect("Loading error");
        assert_eq!(ng.tablesizes(), &[19, 17, 13, 11, 7, 5]);
        assert_eq!(ng.ksize(), 3);
        assert_eq!(ng.get_kmer(b"ACG"), 1);
        assert_eq!(ng.get_kmer(b"TTA"), 1);
        assert_eq!(ng.get_kmer(b"CGA"), 1);
    }

    #[test]
    fn count_and_get_nodegraph() {
        let mut ng: Nodegraph = Nodegraph::new(&[10], 3);

        ng.count(801084876663808);

        assert_eq!(ng.get(801084876663808), 1);
        assert_eq!(ng.unique_kmers(), 1);
    }

    #[test]
    fn load_save_nodegraph() {
        let mut datadir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        datadir.push("../../tests/test-data/.sbt.v3/");

        for i in 0..=5 {
            let mut filename = datadir.clone();
            filename.push(format!("internal.{}", i));
            let data = std::fs::read(filename).unwrap();

            let mut reader = BufReader::new(&data[..]);

            let ng: Nodegraph = Nodegraph::from_reader(&mut reader).expect("Loading error");

            let mut buf = Vec::new();
            {
                let mut writer = BufWriter::new(&mut buf);
                ng.save_to_writer(&mut writer).unwrap();
            }

            let chunk_size = 8;
            for (c1, c2) in data.to_vec().chunks(chunk_size).zip(buf.chunks(chunk_size)) {
                assert_eq!(c1, c2);
            }
            assert_eq!(data.len(), buf.len());
        }
    }

    #[test]
    fn binary_repr_load() {
        let mut reader = BufReader::new(&RAW_DATA[..]);
        let khmer_ng: Nodegraph = Nodegraph::from_reader(&mut reader).expect("Loading error");
        assert_eq!(khmer_ng.tablesizes(), &[19, 17, 13, 11, 7, 5]);
        assert_eq!(khmer_ng.ksize(), 3);
        assert_eq!(khmer_ng.get_kmer(b"ACG"), 1);
        assert_eq!(khmer_ng.get_kmer(b"TTA"), 1);
        assert_eq!(khmer_ng.get_kmer(b"CGA"), 1);

        let mut ng = Nodegraph::with_tables(23, 6, 3);
        ng.count_kmer(b"ACG");
        ng.count_kmer(b"TTA");
        ng.count_kmer(b"CGA");

        let mut buf = Vec::new();
        {
            let mut writer = BufWriter::new(&mut buf);
            ng.save_to_writer(&mut writer).unwrap();
        }
        assert_eq!(buf.len(), 79);
        assert_eq!(&RAW_DATA, &buf.as_slice());
    }

    #[test]
    fn binary_repr_save() {
        let mut ng = Nodegraph::with_tables(23, 6, 3);
        ng.count_kmer(b"ACG");
        ng.count_kmer(b"TTA");
        ng.count_kmer(b"CGA");

        let mut buf = Vec::new();
        {
            let mut writer = BufWriter::new(&mut buf);
            ng.save_to_writer(&mut writer).unwrap();
        }
        let mut reader = BufReader::new(&buf[..]);
        let new_ng: Nodegraph = Nodegraph::from_reader(&mut reader).expect("Loading error");
        assert_eq!(new_ng.tablesizes(), &[19, 17, 13, 11, 7, 5]);
        assert_eq!(new_ng.ksize(), 3);
        assert_eq!(new_ng.get_kmer(b"ACG"), 1);
        assert_eq!(new_ng.get_kmer(b"TTA"), 1);
        assert_eq!(new_ng.get_kmer(b"CGA"), 1);

        assert_eq!(buf.len(), 79);
        assert_eq!(&RAW_DATA, &buf.as_slice());
    }

    #[test]
    fn update_nodegraph() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/.sbt.v3/internal.0");

        let ng_parent: Nodegraph = Nodegraph::from_path(filename).expect("Loading error");

        filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/.sbt.v3/internal.1");

        let ng_1: Nodegraph = Nodegraph::from_path(filename).expect("Loading error");

        filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/.sbt.v3/internal.2");

        let ng_2: Nodegraph = Nodegraph::from_path(filename).expect("Loading error");

        let mut ng_0: Nodegraph = Nodegraph::new(&[99991, 99989, 99971, 99961], 1);
        ng_1.update(&mut ng_0).expect("Error in update");
        ng_2.update(&mut ng_0).expect("Error in update");
        assert_eq!(ng_0.bs, ng_parent.bs);
        //assert_eq!(ng_0.occupied_bins, ng_parent.occupied_bins);
    }

    #[test]
    fn update_nodegraph_many() -> Result<(), Box<dyn std::error::Error>> {
        let mut leaf1 = Nodegraph::with_tables(100, 3, 5);
        for kmer in &["AAAAA", "AAAAT", "AAAAC"] {
            leaf1.count(_hash(kmer.as_bytes()));
        }

        let mut leaf2 = Nodegraph::with_tables(100, 3, 5);
        for kmer in &["AAAAA", "AAAAT", "AAAAG"] {
            leaf2.count(_hash(kmer.as_bytes()));
        }

        let mut leaf3 = Nodegraph::with_tables(100, 3, 5);
        for kmer in &["AAAAA", "AAAAT", "CAAAA"] {
            leaf3.count(_hash(kmer.as_bytes()));
        }

        let mut leaf4 = Nodegraph::with_tables(100, 3, 5);
        for kmer in &["AAAAA", "CAAAA", "GAAAA"] {
            leaf4.count(_hash(kmer.as_bytes()));
        }

        let mut leaf5 = Nodegraph::with_tables(100, 3, 5);
        for kmer in &["AAAAA", "AAAAT", "GAAAA"] {
            leaf5.count(_hash(kmer.as_bytes()));
        }

        let h = _hash(b"AAAAT");
        for leaf in &[leaf1, leaf2, leaf3, leaf5] {
            assert_eq!(leaf.get(h), 1);
        }

        assert_eq!(leaf4.get(h), 0);

        Ok(())
    }

    #[test]
    fn load_nodegraph() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/.sbt.v3/internal.0");

        //let data = include_bytes!("data/internal.0");

        let ng: Nodegraph = Nodegraph::from_path(filename).expect("Loading error");

        assert_eq!(ng.tablesizes(), [99991, 99989, 99971, 99961]);
        assert_eq!(ng.n_occupied_bins(), 2416);
        assert_eq!(ng.get(1877811740), 0);
        for h in [
            1877811749,
            1339603207230,
            5641354835174,
            10502027926594,
            11550845136154,
            12183113567732,
            14062071191653,
            14580861632266,
            18722876140337,
            20714320729467,
            22732389403804,
            24134363957219,
            30606147678309,
            30841792132441,
            31130970675642,
            32760645340554,
            33190965408032,
            33960067474598,
            35413666412010,
            37166860055638,
            38008340488610,
            38631948370393,
            38946626358857,
            39177463395973,
            39396232170068,
            40000457533067,
            41548684950793,
            42975853122398,
            43119393989323,
            43377695911881,
            49367718187361,
            49468277378328,
            50266038601832,
            51636068122286,
            56622962479482,
            58428533496606,
            58971444597606,
            59372670276820,
            59452528403612,
            61074441390615,
            62130354354877,
            62702978264830,
            64430859773984,
            65419869837915,
            65663647257358,
            67872638217057,
            68827108109263,
            69134145403133,
            70436552236751,
            70880519905358,
            78004711377952,
            81502993782978,
            84636365982041,
            85239629151685,
            94266407193778,
            98142256300701,
            98837920540443,
            99930975216128,
            100653760748845,
            102082282949673,
            102530908835648,
            103010972337870,
            103329805967682,
            103652023867250,
            104130252812879,
            112760650992638,
            114779375695317,
            115796389594898,
            117864921668170,
            119763283100790,
            120285237540732,
            121866736124647,
            122140892054804,
            122995254140976,
            123065069359489,
            123405856681590,
            128261346941417,
            130618284885748,
            131310062444107,
            133580282506938,
            139762252968300,
            148434659896290,
            150472163116319,
            151610888790844,
            151736593364935,
            152145317861349,
            154119208822262,
            154803963303860,
            164146490870545,
            166146331478050,
            166719940886532,
            173367021064967,
            173503876669758,
            173949973069402,
            175345218226732,
            175559849681044,
            177057739236298,
            182134979074863,
            185526639726849,
            186188120396587,
            191078441509481,
            191784713609488,
            196150349451960,
            196584209022550,
            196853921592387,
            197752504251580,
            198597053692927,
            200567230796156,
            201179164742411,
            202960515626517,
            203378213499023,
            210822710165852,
            211915017282095,
            213613291536686,
            215418355892998,
            216444054660744,
            216772483699428,
            218586803538885,
            219619606513837,
            221322641419906,
            221692515333150,
            222646058515199,
            223103766020907,
            223436957406949,
            225216425962890,
            225962923363564,
            227026140769845,
            227790244540446,
            228251083676258,
            231710804058239,
            233288106176435,
            235385609463388,
            235438505061770,
            238869764444344,
            239420157045937,
            241121021240187,
            241671335688938,
            242838856557679,
            244786468497109,
            247140303430449,
            248336783901894,
            250357693564448,
            253975323975963,
            256375919657769,
            259301238714261,
            265736169322750,
            265781739304017,
            266725362494513,
            267345873524094,
            271342665825792,
            274876788032658,
            275360996806051,
            275711441656065,
            276221877341287,
            277115529175674,
            277862338800417,
            280967669495427,
            281817613252845,
            281897628539431,
            282200323162036,
            284620358398045,
            284881057128884,
            285925400570356,
            289038917997203,
            289724862541255,
            290309864993733,
            294086384353867,
            295503963521838,
            296966685834878,
            299005107402724,
            300199234365396,
            300617258525997,
            301443933468348,
            302667628736144,
            305781540735975,
            308107503975413,
            308473366560206,
            311148974624393,
            311393227334671,
            312856558437716,
            314634385460120,
            315140251773348,
            316147818305256,
            317314266550052,
            318043998368340,
            319121931997971,
            324333149672473,
            324779561826125,
            326855577904572,
            327646715321140,
            332098363218169,
            333944737799563,
            334160175766170,
            335584394916553,
            335971123608722,
            336472954791992,
            338443948117005,
            338762957149102,
            341091055062112,
            341724341043975,
            343240684449173,
            344010897833199,
            345196014534640,
            347580313704916,
            348815216366639,
            348987115477673,
            350399163507829,
            357535517122796,
            358595265377108,
            358821394913517,
            359452645935849,
            362124977362793,
            366354200059782,
            366535672236781,
            369474755519844,
            370249620342175,
            372037414685096,
            373949557068914,
            374319819178480,
            374609596539290,
            374615513078797,
            375780195152331,
            379102542404949,
            379241504134406,
            379468459802010,
            379661395441316,
            382035531157070,
            383008100523152,
            383135333541903,
            383850900061929,
            384049466048679,
            386263487549463,
            389141313731258,
            390332660259608,
            393516543506060,
            400967959890432,
            401487977714282,
            403579902131163,
            406955472999822,
            408962716867059,
            409903018669983,
            410861197839878,
            414355853800959,
            416580890530128,
            418934773149726,
            419642123579295,
            421963163293847,
            423404494960378,
            424303224424616,
            424596150389604,
            427230335237565,
            429952924284227,
            430664272577516,
            432630098291297,
            434623968464695,
            435267549331128,
            435277763415865,
            435874505125675,
            437654980371254,
            438061138128325,
            438738288109196,
            439177016005977,
            445344075816835,
            445802335759252,
            446710003143163,
            447467518423055,
            449641727299803,
            450058424424520,
            450112320572118,
            450125274173050,
            452241247094714,
            452829154656306,
            454813132622585,
            456174765596578,
            456493632715805,
            456717723773303,
            461156956524045,
            462211497323948,
            463604028403361,
            465228093393002,
            466250095735125,
            469687793491358,
            471922058927200,
            472039595540269,
            472566025949945,
            472595419353109,
            472977022618999,
            473018780652067,
            473772140307174,
            474570287539184,
            474912397870603,
            476325119891604,
            476526896773980,
            476855560317170,
            480232815782455,
            484291524803718,
            485278877010947,
            487732314724511,
            491715999174683,
            494276065129917,
            495846359323641,
            506531113930798,
            507871334392190,
            508031302306958,
            508934816424512,
            509939413858428,
            510737910464301,
            512514768813167,
            513350289212553,
            517460246914282,
            523321188654478,
            524296526109332,
            525762219690878,
            526111205078257,
            527062179866457,
            527591752682839,
            527920198105606,
            530316966667021,
            532977797373940,
            533221992957154,
            533383900955463,
            537527309474265,
            538136383284668,
            538939534540869,
            539777176029418,
            539873986742508,
            543935720187395,
            545273268128445,
            549484636278027,
            551381720133873,
            553977959695484,
            555321949850378,
            555828795847874,
            557285930201258,
            558008777268240,
            558433475619762,
            558892016080993,
            559199414492426,
            560748186311107,
            561604684739024,
            562789967643507,
            563343385252253,
            563775395645616,
            564616206473372,
            565020390122451,
            568901431510366,
            572526115602502,
            573767900523468,
            573851852316852,
            576624529060777,
            576874504697497,
            578856083248351,
            579395263040626,
            579656586099131,
            584217116139474,
            587458649504773,
            591009756408904,
            592792708776319,
            592997432856726,
            594482884410814,
            596004492939074,
            596726606390901,
            597875929908982,
            600179982751750,
            601000534535072,
            601440269988372,
            601603906866038,
            602082770371066,
            604883041984487,
            605545396594434,
            606419362199228,
            607833403537880,
            609555580824872,
            609609500753196,
            611579272742038,
            612206643585093,
            612640334623643,
            612821302220884,
            617021904160724,
            617244669177560,
            617309228629787,
            618709483466270,
            620059729516362,
            620849299055244,
            621083126852990,
            622843084945666,
            623088556560813,
            627738708322473,
            628002002108775,
            628967244202734,
            630034340392901,
            632757066611488,
            634340585739407,
            634691502028135,
            635939425862264,
            637603178700210,
            637880811482435,
            644557275230225,
            644935615624623,
            645793929303122,
            646731502743275,
            646973138978211,
            647900742708077,
            649351154360370,
            653652775436966,
            655230244020599,
            668170744538822,
            670595660720839,
            671785773373187,
            672641554971634,
            672821857332020,
            673587502056476,
            676044446355190,
            677295740685782,
            679716691783353,
            682874745971459,
            682963108550465,
            683897063771844,
            685246440558482,
            686035384279530,
            687129162879229,
            687440351836027,
            688990372747831,
            690608944213791,
            691680901171966,
            694851976547107,
            694869046270466,
            700054088308311,
            701010566680671,
            701156706346414,
            702431887238370,
            702728791577749,
            703127461004015,
            703460523248065,
            705302678110381,
            707793984897058,
            707799855432305,
            707962189637436,
            707993631271976,
            708854130532070,
            710403353214581,
            710927468728191,
            711091480855740,
            712661928452840,
            715334925158742,
            715763419567022,
            715896323316677,
            717568681000032,
            717790011003345,
            719139881875323,
            722537026567926,
            722774506110892,
            723332805980528,
            724621545164802,
            724746920000049,
            727030394121071,
            727262050490847,
            728279662753580,
            730854175545196,
            731361512976697,
            734622692371860,
            736290151677476,
            737921635760471,
            738115824615020,
            739389456325310,
            742704052187442,
            746469097917429,
            748064810280445,
            749144352424687,
            753113822684627,
            753423569783277,
            755196264392026,
            758186007844395,
            758543555642030,
            759083903793759,
            761260029175908,
            767230586289375,
            770167973924874,
            770328708409334,
            772165475523258,
            772947318346532,
            774312511311396,
            774365323868051,
            774964429534347,
            775558532281404,
            779330069525835,
            781344931111517,
            787747218685488,
            788027556261557,
            790211243959626,
            790890494413778,
            792003960897692,
            792629819473398,
            797511060014001,
            797622366845781,
            799257433888961,
            800060479182618,
            801084876663808,
            802340523858506,
            803596407436267,
        ]
        .iter()
        {
            assert_eq!(ng.get(*h), 1);
        }
    }
}
