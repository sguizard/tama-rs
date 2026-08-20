//! `tama split` — file splitting tools.

use clap::{Args as ClapArgs, Subcommand};

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Split a FASTA into chunks. (tama_fasta_splitter)
    Fasta,
    /// Split a mapped SAM into chunks. (tama_mapped_sam_splitter)
    Sam,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let name = match args.cmd {
        Cmd::Fasta => "split fasta",
        Cmd::Sam => "split sam",
    };
    Err(super::not_implemented(name))
}
