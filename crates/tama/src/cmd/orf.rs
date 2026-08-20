//! `tama orf` — ORF / NMD prediction tools.

use clap::{Args as ClapArgs, Subcommand};

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Find ORFs in transcript sequences. (tama_orf_seeker)
    Seek,
    /// Extract CDS regions from a BED. (tama_bed_extract_cds)
    ExtractCds,
    /// Add CDS regions to a BED. (tama_cds_regions_bed_add)
    AddCds,
    /// Parse blastp output for ORF selection. (tama_orf_blastp_parser)
    BlastpParse,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let name = match args.cmd {
        Cmd::Seek => "orf seek",
        Cmd::ExtractCds => "orf extract-cds",
        Cmd::AddCds => "orf add-cds",
        Cmd::BlastpParse => "orf blastp-parse",
    };
    Err(super::not_implemented(name))
}
