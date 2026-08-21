//! Fixed-value CLI options shared by more than one subcommand group.
//!
//! The original TAMA scripts spell these values with underscores (`no_cap`,
//! `common_ends`, …), so every enum here carries
//! `#[value(rename_all = "snake_case")]` — clap's default would render them as
//! `no-cap` / `common-ends` and break compatibility with existing command lines.
//!
//! These are the *CLI* spellings. Where a domain enum already exists in
//! `tama-core`, the conversion lives here so `tama-core` stays free of any clap
//! dependency.

use clap::ValueEnum;

use tama_core::collapse::{Cap, Ends};
use tama_core::metrics::IdentMethod;

/// 5' cap handling. (`-x`)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[value(rename_all = "snake_case")]
pub enum CapFlag {
    /// 5' cap-trimmed data; 5' ends are trusted.
    Capped,
    /// Non-cap-trimmed data; 5' degradation is expected.
    NoCap,
}

impl From<CapFlag> for Cap {
    fn from(v: CapFlag) -> Self {
        match v {
            CapFlag::Capped => Cap::Capped,
            CapFlag::NoCap => Cap::NoCap,
        }
    }
}

/// Exon-end collapsing behaviour. (`-e`)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[value(rename_all = "snake_case")]
pub enum EndsOpt {
    /// Use the most common start/end among the collapsed models.
    CommonEnds,
    /// Use the longest start/end among the collapsed models.
    LongestEnds,
}

impl From<EndsOpt> for Ends {
    fn from(v: EndsOpt) -> Self {
        match v {
            EndsOpt::CommonEnds => Ends::CommonEnds,
            EndsOpt::LongestEnds => Ends::LongestEnds,
        }
    }
}

/// Identity calculation method. (`-icm`)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[value(rename_all = "snake_case")]
pub enum IdentMethodOpt {
    /// Include hard/soft clipping in the denominator (default).
    IdentCov,
    /// Exclude hard/soft clipping.
    IdentMap,
}

impl From<IdentMethodOpt> for IdentMethod {
    fn from(v: IdentMethodOpt) -> Self {
        match v {
            IdentMethodOpt::IdentCov => IdentMethod::IdentCov,
            IdentMethodOpt::IdentMap => IdentMethod::IdentMap,
        }
    }
}

/// Duplicate-model merge behaviour. (`-d`)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[value(rename_all = "snake_case")]
pub enum Dup {
    /// Merge duplicate models into one.
    MergeDup,
    /// Keep duplicate models separate.
    NoMerge,
}

/// Splice-junction priority. (`-sj`)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[value(rename_all = "snake_case")]
pub enum SjPriority {
    /// Do not prioritise splice junctions over ends.
    NoPriority,
    /// Prioritise splice junctions over ends.
    SjPriority,
}

/// Collapse run mode. (`-rm`)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[value(rename_all = "snake_case")]
pub enum RunMode {
    /// Hold everything in memory (default).
    Original,
    /// Lower-memory mode.
    LowMem,
}

/// Whether filtering removes whole genes or single transcripts. (`-l`)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[value(rename_all = "snake_case")]
pub enum Level {
    /// Remove at the gene level.
    Gene,
    /// Remove at the transcript level.
    Transcript,
}

/// Multi-exon model handling during filtering. (`-k`)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[value(rename_all = "snake_case")]
pub enum Multi {
    /// Keep multi-exon models regardless of support.
    KeepMulti,
    /// Subject multi-exon models to the same filter.
    RemoveMulti,
}

/// Whether the stop codon is part of the CDS. (`-s`)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[value(rename_all = "snake_case")]
pub enum StopCodon {
    /// Include the stop codon in the CDS.
    IncludeStop,
    /// Exclude the stop codon from the CDS.
    NoStopCodon,
}
