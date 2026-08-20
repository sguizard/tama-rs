//! `tama stats` — file statistics tools.

use std::io::Write;

use clap::{Args as ClapArgs, Subcommand};
use indexmap::IndexMap;

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
    /// Sampling saturation curve from a read-support levels file. (tama_sampling_saturation_curve)
    Saturation {
        /// Read support levels file. (`-r`)
        #[arg(short = 'r', long)]
        report: std::path::PathBuf,
        /// Read bin size. (`-b`)
        #[arg(short = 'b', long)]
        bin: usize,
        /// Output file. (`-o`)
        #[arg(short = 'o', long)]
        output: std::path::PathBuf,
    },
}

pub fn run(args: Args) -> anyhow::Result<()> {
    match args.cmd {
        Cmd::Degradation => Err(super::not_implemented("stats degradation")),
        Cmd::ModelChanges => Err(super::not_implemented("stats model-changes")),
        Cmd::Saturation { report, bin, output } => saturation(&report, bin, &output),
    }
}

/// Sampling saturation curve. Ports `tama_sampling_saturation_curve`.
///
/// Note: the curve is cumulative over reads in the order they first appear; the
/// original iterated a Python-2 dict (hash order), so only the final totals are
/// guaranteed to match — the intermediate sampling points differ.
fn saturation(report: &std::path::Path, bin: usize, output: &std::path::Path) -> anyhow::Result<()> {
    use std::io::BufRead;
    let reader = tama_io::open_reader(report)?;
    let mut read_gene: IndexMap<String, String> = IndexMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with("gene_id") || line.starts_with("merge_gene_id") {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 6 {
            continue;
        }
        let gene_id = cols[0];
        for source_line in cols[5].split(';') {
            let parts: Vec<&str> = source_line.split(':').collect();
            if parts.len() < 2 {
                continue;
            }
            let read_line = parts[1..].join(":");
            for read_id in read_line.split(',') {
                read_gene
                    .entry(read_id.to_string())
                    .or_insert_with(|| gene_id.to_string());
            }
        }
    }

    let mut out = tama_io::create_writer(output)?;
    writeln!(out, "read_count\tgene_count")?;
    let mut read_count = 0usize;
    let mut gene_count = 0usize;
    let mut seen: std::collections::HashSet<&str> = std::collections::HashSet::new();
    for (read_id, gene) in &read_gene {
        read_count += 1;
        if bin != 0 && read_count % bin == 0 {
            writeln!(out, "{read_count}\t{gene_count}")?;
        }
        let _ = read_id;
        if seen.insert(gene.as_str()) {
            gene_count += 1;
        }
    }
    Ok(())
}
