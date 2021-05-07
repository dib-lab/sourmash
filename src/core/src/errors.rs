use thiserror::Error;

#[derive(Debug, Error)]
pub enum SourmashError {
    /// Raised for internal errors in the libraries.  Should not happen.
    #[error("internal error: {message:?}")]
    Internal { message: String },

    #[error("must have same num: {n1} != {n2}")]
    MismatchNum { n1: u32, n2: u32 },

    #[error("different ksizes cannot be compared")]
    MismatchKSizes,

    #[error("DNA/prot minhashes cannot be compared")]
    MismatchDNAProt,

    #[error("mismatch in scaled; comparison fail")]
    MismatchScaled,

    #[error("mismatch in seed; comparison fail")]
    MismatchSeed,

    #[error("different signatures cannot be compared")]
    MismatchSignatureType,

    #[error("Invalid hash function: {function:?}")]
    InvalidHashFunction { function: String },

    #[error("Can only set {message:?} if the MinHash is empty")]
    NonEmptyMinHash { message: String },

    #[error("invalid DNA character in input k-mer: {message}")]
    InvalidDNA { message: String },

    #[error("invalid protein character in input: {message}")]
    InvalidProt { message: String },

    #[error("Codon is invalid length: {message}")]
    InvalidCodonLength { message: String },

    #[error("Set error rate to a value smaller than 0.367696 and larger than 0.00203125")]
    HLLPrecisionBounds,

    #[error(transparent)]
    ReadDataError(#[from] crate::storage::ReadDataError),

    #[error(transparent)]
    StorageError(#[from] crate::storage::StorageError),

    #[error(transparent)]
    SerdeError(#[from] serde_json::error::Error),

    #[error(transparent)]
    NifflerError(#[from] niffler::Error),

    #[error(transparent)]
    Utf8Error(#[from] std::str::Utf8Error),

    #[error(transparent)]
    IOError(#[from] std::io::Error),

    #[cfg(not(all(target_arch = "wasm32", target_vendor = "unknown")))]
    #[error(transparent)]
    Panic(#[from] crate::ffi::utils::Panic),
}

#[repr(u32)]
pub enum SourmashErrorCode {
    // no error
    NoError = 0,
    // panics and internals
    Panic = 1,
    Internal = 2,
    Msg = 3,
    Unknown = 4,
    // Compatibility errors
    MismatchKSizes = 1_01,
    MismatchDNAProt = 1_02,
    MismatchScaled = 1_03,
    MismatchSeed = 1_04,
    MismatchSignatureType = 1_05,
    NonEmptyMinHash = 1_06,
    MismatchNum = 1_07,
    // Input sequence errors
    InvalidDNA = 11_01,
    InvalidProt = 11_02,
    InvalidCodonLength = 11_03,
    InvalidHashFunction = 11_04,
    // index-related errors
    ReadData = 12_01,
    Storage = 12_02,
    // HLL errors
    HLLPrecisionBounds = 13_01,
    // external errors
    Io = 100_001,
    Utf8Error = 100_002,
    ParseInt = 100_003,
    SerdeError = 100_004,
    NifflerError = 100_005,
}

#[cfg(not(all(target_arch = "wasm32", target_vendor = "unknown")))]
impl SourmashErrorCode {
    pub fn from_error(error: &SourmashError) -> SourmashErrorCode {
        match error {
            SourmashError::Internal { .. } => SourmashErrorCode::Internal,
            SourmashError::Panic { .. } => SourmashErrorCode::Panic,
            SourmashError::MismatchNum { .. } => SourmashErrorCode::MismatchNum,
            SourmashError::MismatchKSizes => SourmashErrorCode::MismatchKSizes,
            SourmashError::MismatchDNAProt => SourmashErrorCode::MismatchDNAProt,
            SourmashError::MismatchScaled => SourmashErrorCode::MismatchScaled,
            SourmashError::MismatchSeed => SourmashErrorCode::MismatchSeed,
            SourmashError::MismatchSignatureType => SourmashErrorCode::MismatchSignatureType,
            SourmashError::NonEmptyMinHash { .. } => SourmashErrorCode::NonEmptyMinHash,
            SourmashError::InvalidDNA { .. } => SourmashErrorCode::InvalidDNA,
            SourmashError::InvalidProt { .. } => SourmashErrorCode::InvalidProt,
            SourmashError::InvalidCodonLength { .. } => SourmashErrorCode::InvalidCodonLength,
            SourmashError::InvalidHashFunction { .. } => SourmashErrorCode::InvalidHashFunction,
            SourmashError::ReadDataError { .. } => SourmashErrorCode::ReadData,
            SourmashError::StorageError { .. } => SourmashErrorCode::Storage,
            SourmashError::HLLPrecisionBounds { .. } => SourmashErrorCode::HLLPrecisionBounds,
            SourmashError::SerdeError { .. } => SourmashErrorCode::SerdeError,
            SourmashError::IOError { .. } => SourmashErrorCode::Io,
            SourmashError::NifflerError { .. } => SourmashErrorCode::NifflerError,
            SourmashError::Utf8Error { .. } => SourmashErrorCode::Utf8Error,
        }
    }
}
