pub mod hyperloglog;
pub mod minhash;

#[cfg(not(all(target_arch = "wasm32", target_vendor = "unknown")))]
pub mod nodegraph;

use serde::{Deserialize, Serialize};

use crate::sketch::hyperloglog::HyperLogLog;
use crate::sketch::minhash::{KmerMinHash, KmerMinHashBTree};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Sketch {
    MinHash(KmerMinHash),
    LargeMinHash(KmerMinHashBTree),
    HyperLogLog(HyperLogLog),
}
