//! `tama stats` — file statistics tools.

use clap::{Args as ClapArgs, Subcommand};

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Degradation signature. (tama_degradation_signature)
    Degradation,
    /// Find model changes between annotations. (tama_find_model_changes)
    ModelChanges,
    /// Sampling saturation curve. (tama_sampling_saturation_curve)
    Saturation,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let name = match args.cmd {
        Cmd::Degradation => "stats degradation",
        Cmd::ModelChanges => "stats model-changes",
        Cmd::Saturation => "stats saturation",
    };
    Err(super::not_implemented(name))
}
