//! `tama` — Transcriptome Annotation by Modular Algorithms (Rust rewrite).
//!
//! A single binary exposing every original TAMA tool as a subcommand. Run
//! `tama <group> --help` to see the tools in each group.

use clap::{Parser, Subcommand};
use tama::cmd;

#[derive(Parser)]
#[command(
    name = "tama",
    version,
    about = "Transcriptome Annotation by Modular Algorithms",
    long_about = "TAMA processes long-read (Iso-Seq / Nanopore) transcriptome data. \
                  Each original TAMA script is available here as a subcommand."
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Collapse mapped reads into transcript models (tama_collapse).
    Collapse(cmd::collapse::Args),
    /// Merge transcript annotations across sources (tama_merge).
    Merge(cmd::merge::Args),
    /// Format converters (BED <-> GTF/GFF, FASTQ -> FASTA, ...).
    Format(cmd::format::Args),
    /// Filter transcript models.
    Filter(cmd::filter::Args),
    /// ORF / NMD prediction tools.
    Orf(cmd::orf::Args),
    /// Read-support tracking tools.
    Support(cmd::support::Args),
    /// File statistics tools.
    Stats(cmd::stats::Args),
    /// File splitting tools.
    Split(cmd::split::Args),
    /// Variant calling.
    Variants(cmd::variants::Args),
    /// Sequence cleanup tools.
    Cleanup(cmd::cleanup::Args),
}

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp(None)
        .init();

    let cli = Cli::parse();
    match cli.command {
        Command::Collapse(a) => cmd::collapse::run(a),
        Command::Merge(a) => cmd::merge::run(a),
        Command::Format(a) => cmd::format::run(a),
        Command::Filter(a) => cmd::filter::run(a),
        Command::Orf(a) => cmd::orf::run(a),
        Command::Support(a) => cmd::support::run(a),
        Command::Stats(a) => cmd::stats::run(a),
        Command::Split(a) => cmd::split::run(a),
        Command::Variants(a) => cmd::variants::run(a),
        Command::Cleanup(a) => cmd::cleanup::run(a),
    }
}
