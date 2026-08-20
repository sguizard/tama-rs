//! `tama collapse` — collapse mapped long reads into transcript models.
//!
//! Ports `tama_collapse.py`. Argument names mirror the original short flags so
//! existing pipelines translate directly.

use clap::Parser;

#[derive(Parser)]
pub struct Args {
    /// Sorted SAM file (or BAM with --bam). (`-s`)
    #[arg(short = 's', long = "sam")]
    pub sam: std::path::PathBuf,

    /// Genome FASTA file. (`-f`)
    #[arg(short = 'f', long = "fasta")]
    pub fasta: std::path::PathBuf,

    /// Output prefix. (`-p`)
    #[arg(short = 'p', long = "prefix")]
    pub prefix: String,

    /// Capped flag: `capped` or `no_cap`. (`-x`)
    #[arg(short = 'x', long, default_value = "no_cap")]
    pub cap_flag: String,

    /// Collapse exon ends: `common_ends` or `longest_ends`. (`-e`)
    #[arg(short = 'e', long, default_value = "common_ends")]
    pub ends: String,

    /// Minimum coverage percent. (`-c`)
    #[arg(short = 'c', long, default_value_t = 99.0)]
    pub coverage: f64,

    /// Minimum identity percent. (`-i`)
    #[arg(short = 'i', long, default_value_t = 85.0)]
    pub identity: f64,

    /// Identity calculation method: `ident_cov` or `ident_map`. (`-icm`)
    #[arg(long = "icm", default_value = "ident_cov")]
    pub ident_method: String,

    /// 5' threshold. (`-a`)
    #[arg(short = 'a', long, default_value_t = 10)]
    pub five_prime: i64,

    /// Exon/splice-junction threshold. (`-m`)
    #[arg(short = 'm', long, default_value_t = 10)]
    pub exon_thresh: i64,

    /// 3' threshold. (`-z`)
    #[arg(short = 'z', long, default_value_t = 10)]
    pub three_prime: i64,

    /// Duplicate merge behaviour: `merge_dup` or `no_merge`. (`-d`)
    #[arg(short = 'd', long, default_value = "merge_dup")]
    pub dup: String,

    /// Splice-junction priority: `no_priority` or `sj_priority`. (`-sj`)
    #[arg(long = "sj", default_value = "no_priority")]
    pub sj_priority: String,

    /// Splice-junction error threshold (bp). (`-sjt`)
    #[arg(long = "sjt", default_value_t = 10)]
    pub sj_thresh: i64,

    /// Local density error threshold. (`-lde`)
    #[arg(long = "lde", default_value_t = 1000)]
    pub lde: i64,

    /// Simple error symbol for LDE output. (`-ses`)
    #[arg(long = "ses", default_value = "^")]
    pub simple_error_symbol: String,

    /// Treat input as BAM instead of SAM. (`-b`)
    #[arg(short = 'b', long)]
    pub bam: bool,

    /// Run mode: `original` or `low_mem`. (`-rm`)
    #[arg(long = "rm", default_value = "original")]
    pub run_mode: String,

    /// Variation coverage threshold (reads). (`-vc`)
    #[arg(long = "vc", default_value_t = 5)]
    pub var_coverage: i64,
}

pub fn run(_args: Args) -> anyhow::Result<()> {
    Err(super::not_implemented("tama collapse"))
}
