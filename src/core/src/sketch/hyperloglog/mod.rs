/*
Based on the HyperLogLog implementations in khmer
  https://github.com/dib-lab/khmer/blob/fb65d21eaedf0d397d49ae3debc578897f9d6eb4/src/oxli/hllcounter.cc
using the maximum likelihood estimators from
  https://oertl.github.io/hyperloglog-sketch-estimation-paper/paper/paper.pdf
first implemented for genomics in dashing
  https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1875-0
*/

use std::cmp;
use std::fs::File;
use std::io;
use std::path::Path;

use byteorder::{BigEndian, ReadBytesExt, WriteBytesExt};
use serde::{Deserialize, Serialize};

use crate::encodings::HashFunctions;
use crate::signature::SigsTrait;
use crate::sketch::KmerMinHash;
use crate::traits::Update;
use crate::Error;
use crate::HashIntoType;

pub mod estimators;
use estimators::CounterType;

#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
pub struct HyperLogLog {
    registers: Vec<CounterType>,
    p: usize,
    q: usize,
    ksize: usize,
}

impl HyperLogLog {
    pub fn with_error_rate(error_rate: f64, ksize: usize) -> Result<HyperLogLog, Error> {
        let p = f64::ceil(f64::log2(f64::powi(1.04 / error_rate, 2)));
        HyperLogLog::new(p as usize, ksize)
    }

    pub fn new(p: usize, ksize: usize) -> Result<HyperLogLog, Error> {
        if !(4..=18).contains(&p) {
            return Err(Error::HLLPrecisionBounds);
        }

        let size = (1_usize) << p;
        let registers = vec![0; size];

        Ok(HyperLogLog {
            registers,
            ksize,
            p,
            q: 64 - p, // FIXME: allow setting q explicitly
        })
    }

    pub fn merge(&mut self, other: &HyperLogLog) -> Result<(), Error> {
        self.check_compatible(other)?;
        self.registers
            .iter_mut()
            .zip(other.registers.iter())
            .for_each(|(a, b)| *a = cmp::max(*a, *b));
        Ok(())
    }

    pub fn add_word(&mut self, word: &[u8]) {
        let hash = crate::_hash_murmur(word, 42); // TODO: decide on seed
        self.add_hash(hash);
    }

    pub fn add_many(&mut self, hashes: &[HashIntoType]) -> Result<(), Error> {
        for min in hashes {
            self.add_hash(*min);
        }
        Ok(())
    }

    pub fn cardinality(&self) -> usize {
        let counts = estimators::counts(&self.registers, self.q);

        estimators::mle(&counts, self.p, self.q, 0.01) as usize
    }

    pub fn similarity(&self, other: &HyperLogLog) -> f64 {
        let (only_a, only_b, intersection) =
            estimators::joint_mle(&self.registers, &other.registers, self.p, self.q);

        intersection as f64 / (only_a + only_b + intersection) as f64
    }

    pub fn containment(&self, other: &HyperLogLog) -> f64 {
        let (only_a, _, intersection) =
            estimators::joint_mle(&self.registers, &other.registers, self.p, self.q);

        intersection as f64 / (only_a + intersection) as f64
    }

    pub fn intersection(&self, other: &HyperLogLog) -> usize {
        let (_, _, intersection) =
            estimators::joint_mle(&self.registers, &other.registers, self.p, self.q);

        intersection
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
        wtr.write_all(b"HLL")?;
        wtr.write_u8(1)?; // version
        wtr.write_u8(self.p as u8)?; // number of bits used for indexing
        wtr.write_u8(self.q as u8)?; // number of bits used for counting leading zeroes
        wtr.write_u8(self.ksize as u8)?; // ksize
        wtr.write_all(&self.registers.as_slice())?;

        Ok(())
    }

    pub fn from_reader<R>(rdr: R) -> Result<HyperLogLog, Error>
    where
        R: io::Read,
    {
        let (mut rdr, _format) = niffler::get_reader(Box::new(rdr))?;

        let signature = rdr.read_u24::<BigEndian>()?;
        assert_eq!(signature, 0x484c4c);

        let version = rdr.read_u8()?;
        assert_eq!(version, 1);

        let p = rdr.read_u8()? as usize;
        let q = rdr.read_u8()? as usize;

        let ksize = rdr.read_u8()? as usize;
        let n_registers = 1 << p;

        let mut registers = vec![0u8; n_registers];
        rdr.read_exact(&mut registers)?;

        Ok(HyperLogLog {
            p,
            q,
            ksize,
            registers,
        })
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<HyperLogLog, Error> {
        let mut reader = io::BufReader::new(File::open(path)?);
        HyperLogLog::from_reader(&mut reader)
    }
}

impl SigsTrait for HyperLogLog {
    fn size(&self) -> usize {
        self.registers.len()
    }

    fn to_vec(&self) -> Vec<u64> {
        self.registers.iter().map(|x| *x as u64).collect()
    }

    fn ksize(&self) -> usize {
        self.ksize as usize
    }

    fn seed(&self) -> u64 {
        // TODO: support other seeds
        42
    }

    fn hash_function(&self) -> HashFunctions {
        //TODO support other hash functions
        HashFunctions::murmur64_DNA
    }

    fn add_hash(&mut self, hash: HashIntoType) {
        let value = hash >> self.p;
        let index = (hash - (value << self.p)) as usize;

        let leftmost = value.leading_zeros() + 1 - (self.p as u32);

        let old_value = self.registers[index];
        self.registers[index] = cmp::max(old_value, leftmost as CounterType);
    }

    fn check_compatible(&self, other: &HyperLogLog) -> Result<(), Error> {
        if self.ksize() != other.ksize() {
            Err(Error::MismatchKSizes)
        } else if self.size() != other.size() {
            // TODO: create new error
            Err(Error::MismatchNum {
                n1: self.size() as u32,
                n2: other.size() as u32,
            })
        } else {
            Ok(())
        }
    }
}

impl Update<HyperLogLog> for KmerMinHash {
    fn update(&self, other: &mut HyperLogLog) -> Result<(), Error> {
        for h in self.mins() {
            other.add_hash(h);
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashSet;
    use std::io::{BufReader, BufWriter, Read};
    use std::path::PathBuf;

    use crate::signature::SigsTrait;
    use needletail::{parse_fastx_file, parse_fastx_reader, Sequence};

    use super::HyperLogLog;

    // TODO: pull more tests from khmer HLL

    #[test]
    fn hll_add() {
        const ERR_RATE: f64 = 0.01;
        const N_UNIQUE: usize = 3356;
        const KSIZE: u8 = 21;

        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/ecoli.genes.fna");

        let mut hll = HyperLogLog::with_error_rate(ERR_RATE, KSIZE as usize).unwrap();
        let mut counter: HashSet<Vec<u8>> = HashSet::new();

        let mut parser = parse_fastx_file(filename).unwrap();
        while let Some(record) = parser.next() {
            let record = record.unwrap();
            let norm_seq = record.normalize(false);
            let rc = norm_seq.reverse_complement();

            hll.add_sequence(&norm_seq, false).unwrap();
            for (_, kmer, _) in norm_seq.canonical_kmers(KSIZE, &rc) {
                counter.insert(kmer.into());
            }
        }

        assert_eq!(counter.len(), N_UNIQUE);

        let abs_error = (1. - (hll.cardinality() as f64 / N_UNIQUE as f64)).abs();
        assert!(abs_error < ERR_RATE, "{}", abs_error);
    }

    #[test]
    fn hll_joint_mle() {
        const ERR_RATE: f64 = 0.01;
        const KSIZE: u8 = 21;

        const N_UNIQUE_H1: usize = 500741;
        const N_UNIQUE_H2: usize = 995845;
        const N_UNIQUE_U: usize = 995845;

        const SIMILARITY: f64 = 0.502783;
        const CONTAINMENT_H1: f64 = 1.;
        const CONTAINMENT_H2: f64 = 0.502783;

        const INTERSECTION: usize = 500838;

        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/genome-s10.fa.gz");

        let mut hll1 = HyperLogLog::with_error_rate(ERR_RATE, KSIZE as usize).unwrap();
        let mut hll2 = HyperLogLog::with_error_rate(ERR_RATE, KSIZE as usize).unwrap();
        let mut hllu = HyperLogLog::with_error_rate(ERR_RATE, KSIZE as usize).unwrap();

        let mut buf = vec![];
        let (mut reader, _) = niffler::from_path(filename).unwrap();
        reader.read_to_end(&mut buf).unwrap();

        let mut parser = parse_fastx_reader(&buf[..]).unwrap();
        while let Some(record) = parser.next() {
            let record = record.unwrap();
            let norm_seq = record.normalize(false);

            hll1.add_sequence(&norm_seq, false).unwrap();
            hllu.add_sequence(&norm_seq, false).unwrap();
        }

        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/genome-s10+s11.fa.gz");

        let mut buf = vec![];
        let (mut reader, _) = niffler::from_path(filename).unwrap();
        reader.read_to_end(&mut buf).unwrap();

        let mut parser = parse_fastx_reader(&buf[..]).unwrap();
        while let Some(record) = parser.next() {
            let record = record.unwrap();
            let norm_seq = record.normalize(false);

            hll2.add_sequence(&norm_seq, false).unwrap();
            hllu.add_sequence(&norm_seq, false).unwrap();
        }

        let abs_error = (1. - (hll1.cardinality() as f64 / N_UNIQUE_H1 as f64)).abs();
        assert!(abs_error < ERR_RATE, "{}", abs_error);

        let abs_error = (1. - (hll2.cardinality() as f64 / N_UNIQUE_H2 as f64)).abs();
        assert!(abs_error < ERR_RATE, "{}", abs_error);

        let similarity = hll1.similarity(&hll2);
        let abs_error = (1. - (similarity / SIMILARITY as f64)).abs();
        assert!(abs_error < ERR_RATE, "{} {}", similarity, SIMILARITY);

        let containment = hll1.containment(&hll2);
        let abs_error = (1. - (containment / CONTAINMENT_H1 as f64)).abs();
        assert!(abs_error < ERR_RATE, "{} {}", containment, CONTAINMENT_H1);

        let containment = hll2.containment(&hll1);
        let abs_error = (1. - (containment / CONTAINMENT_H2 as f64)).abs();
        assert!(abs_error < ERR_RATE, "{} {}", containment, CONTAINMENT_H2);

        let intersection = hll1.intersection(&hll2) as f64;
        let abs_error = (1. - (intersection / INTERSECTION as f64)).abs();
        assert!(abs_error < ERR_RATE, "{} {}", intersection, INTERSECTION);

        hll1.merge(&hll2).unwrap();

        let abs_error = (1. - (hllu.similarity(&hll1) as f64 / 1.)).abs();
        assert!(abs_error < ERR_RATE, "{}", abs_error);

        let abs_error = (1. - (hllu.containment(&hll1) as f64 / 1.)).abs();
        assert!(abs_error < ERR_RATE, "{}", abs_error);

        let abs_error = (1. - (hll1.containment(&hllu) as f64 / 1.)).abs();
        assert!(abs_error < ERR_RATE, "{}", abs_error);

        let intersection = hll1.intersection(&hllu) as f64;
        let abs_error = (1. - (intersection / N_UNIQUE_U as f64)).abs();
        assert!(abs_error < ERR_RATE, "{} {}", intersection, N_UNIQUE_U);
    }

    #[test]
    fn save_load_hll() {
        let mut hll = HyperLogLog::with_error_rate(0.01, 1).expect("error building HLL");
        for i in 1..5000 {
            hll.add_hash(i)
        }

        let mut buf = Vec::new();
        {
            let mut writer = BufWriter::new(&mut buf);
            hll.save_to_writer(&mut writer).unwrap();
        }

        let mut reader = BufReader::new(&buf[..]);
        let hll_new: HyperLogLog = HyperLogLog::from_reader(&mut reader).expect("Loading error");

        assert_eq!(hll_new.p, hll.p);
        assert_eq!(hll_new.q, hll.q);
        assert_eq!(hll_new.registers, hll.registers);
        assert_eq!(hll_new.ksize, hll.ksize);
    }
}
