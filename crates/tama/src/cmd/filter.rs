//! `tama filter` — transcript model filters.

use clap::{Args as ClapArgs, Subcommand};

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Keep primary transcripts by ORF. (tama_filter_primary_transcripts_orf)
    PrimaryOrf,
    /// Remove fragment models. (tama_remove_fragment_models)
    Fragments,
    /// Remove poly-A models by level. (tama_remove_polya_models_levels)
    Polya,
    /// Remove single-read models by level. (tama_remove_single_read_models_levels)
    SingleRead,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let name = match args.cmd {
        Cmd::PrimaryOrf => "filter primary-orf",
        Cmd::Fragments => "filter fragments",
        Cmd::Polya => "filter polya",
        Cmd::SingleRead => "filter single-read",
    };
    Err(super::not_implemented(name))
}
