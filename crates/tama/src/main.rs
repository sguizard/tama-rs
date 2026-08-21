//! `tama` — Transcriptome Annotation by Modular Algorithms (Rust rewrite).
//!
//! A single binary exposing every original TAMA tool as a subcommand. Run
//! `tama <group> --help` to see the tools in each group.

use std::io::Write;

use clap::{CommandFactory, Parser, Subcommand, ValueEnum};
use clap_complete::{generate, Shell};
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
    /// Print progress messages to stderr while running (so you can see the run is
    /// still going). Off by default; only end-of-run summaries are shown. Can also
    /// be controlled with the RUST_LOG environment variable.
    #[arg(short = 'v', long, global = true)]
    verbose: bool,
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
    /// Generate a shell completion script on stdout.
    ///
    /// Write it to your shell's completions directory, e.g.
    /// `tama completions fish > ~/.config/fish/completions/tama.fish`.
    Completions {
        /// Shell to generate completions for.
        #[arg(value_enum)]
        shell: Shell,
    },
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    // Default: `info` (end-of-run summaries only). `--verbose` raises it to `debug`
    // so the per-step progress heartbeats show. RUST_LOG still overrides both.
    let default_level = if cli.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(default_level))
        .format_timestamp(None)
        .init();

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
        Command::Completions { shell } => {
            let mut cmd = Cli::command();
            let name = cmd.get_name().to_string();
            let mut out = std::io::stdout();
            generate(shell, &mut cmd, name.clone(), &mut out);
            if shell == Shell::Fish {
                // clap_complete's fish generator only emits completions for
                // options, not for positional arguments — without this,
                // `tama completions <TAB>` falls back to listing files.
                let shells: Vec<String> = Shell::value_variants()
                    .iter()
                    .map(|s| s.to_string())
                    .collect();
                writeln!(
                    out,
                    "complete -c {name} -n \"__fish_{name}_using_subcommand completions\" -f -a \"{}\"",
                    shells.join(" ")
                )?;
            }
            Ok(())
        }
    }
}
