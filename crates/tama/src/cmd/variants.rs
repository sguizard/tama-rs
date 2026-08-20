//! `tama variants` — variant calling.

use clap::{Args as ClapArgs, Subcommand};

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Call variants from collapse output. (tama_variant_caller)
    Call,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let name = match args.cmd {
        Cmd::Call => "variants call",
    };
    Err(super::not_implemented(name))
}
