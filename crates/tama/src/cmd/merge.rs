//! `tama merge` — merge transcript annotations across sources.
//!
//! Ports `tama_merge.py`.

use clap::Parser;

#[derive(Parser)]
pub struct Args {
    /// File list describing the annotations to merge. (`-f`)
    #[arg(short = 'f', long = "filelist")]
    pub filelist: std::path::PathBuf,

    /// Output prefix. (`-p`)
    #[arg(short = 'p', long = "prefix")]
    pub prefix: String,

    /// Collapse exon ends: `common_ends` or `longest_ends`. (`-e`)
    #[arg(short = 'e', long, default_value = "common_ends")]
    pub ends: String,

    /// 5' threshold. (`-a`)
    #[arg(short = 'a', long, default_value_t = 10)]
    pub five_prime: i64,

    /// Exon/splice-junction threshold. (`-m`)
    #[arg(short = 'm', long, default_value_t = 10)]
    pub exon_thresh: i64,

    /// 3' threshold. (`-z`)
    #[arg(short = 'z', long, default_value_t = 10)]
    pub three_prime: i64,

    /// Duplicate merge behaviour: `no_merge` or `merge_dup`. (`-d`)
    #[arg(short = 'd', long, default_value = "no_merge")]
    pub dup: String,

    /// Use gene/transcript IDs from this merge source. (`-s`)
    #[arg(short = 's', long)]
    pub source_id: Option<String>,

    /// Use CDS from this merge source. (`-cds`)
    #[arg(long = "cds")]
    pub cds_source: Option<String>,
}

pub fn run(_args: Args) -> anyhow::Result<()> {
    Err(super::not_implemented("tama merge"))
}
