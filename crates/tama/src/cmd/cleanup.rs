//! `tama cleanup` — sequence cleanup tools.

use clap::{Args as ClapArgs, Subcommand};

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Clean up poly-A tails from FLNC reads. (tama_flnc_polya_cleanup)
    Polya,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let name = match args.cmd {
        Cmd::Polya => "cleanup polya",
    };
    Err(super::not_implemented(name))
}
