//! Core domain model and algorithms for the TAMA toolkit.
//!
//! This crate holds the format-agnostic pieces shared by every subcommand: the
//! [`model::BedTranscript`] type and TAMA BED12 dialect, CIGAR handling, and
//! sequence utilities. Higher-level algorithms (collapse/merge comparison,
//! error/variation calling) are layered on top in later phases.

pub mod bed;
pub mod cigar;
pub mod error;
pub mod error_calc;
pub mod gene;
pub mod metrics;
pub mod model;
pub mod polya;
pub mod seq;

pub use error::Error;
pub use model::BedTranscript;
