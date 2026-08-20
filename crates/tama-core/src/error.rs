//! Error type shared across the TAMA core algorithms.

use thiserror::Error as ThisError;

#[derive(ThisError, Debug)]
pub enum Error {
    #[error("BED parse error: {0}")]
    Bed(String),

    #[error("CIGAR parse error: {0}")]
    Cigar(String),

    #[error("invalid input: {0}")]
    Invalid(String),

    #[error(transparent)]
    Io(#[from] std::io::Error),
}
