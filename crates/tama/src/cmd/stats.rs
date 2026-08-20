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
    /// Degradation signature from capped vs no_cap collapse. (tama_degradation_signature)
    Degradation {
        /// Capped collapse `trans_read.bed`. (`-c`)
        #[arg(short = 'c', long)]
        capped: std::path::PathBuf,
        /// No-cap collapse `trans_read.bed`. (`-nc`)
        #[arg(long = "nc")]
        nocap: std::path::PathBuf,
        /// Output file. (`-o`)
        #[arg(short = 'o', long)]
        output: std::path::PathBuf,
    },
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
        Cmd::Degradation { capped, nocap, output } => degradation(&capped, &nocap, &output),
        Cmd::ModelChanges => Err(super::not_implemented("stats model-changes")),
        Cmd::Saturation { report, bin, output } => saturation(&report, bin, &output),
    }
}

/// Per-source gene/transcript exon and read tallies.
#[derive(Default)]
struct DegStats {
    gene_order: Vec<String>,
    gene_max_exons: IndexMap<String, usize>,
    gene_reads: IndexMap<String, std::collections::HashSet<String>>,
    gene_trans: IndexMap<String, std::collections::HashSet<String>>,
    trans_max_exons: IndexMap<String, usize>,
}

fn read_deg(path: &std::path::Path) -> anyhow::Result<DegStats> {
    use std::io::BufRead;
    let mut s = DegStats::default();
    for line in tama_io::open_reader(path)?.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        let id_split: Vec<&str> = cols[3].split(';').collect();
        let trans_id = id_split[0].to_string();
        let read_id = id_split[1].to_string();
        let gene_id = trans_id.split('.').next().unwrap_or(&trans_id).to_string();
        let num_exons: usize = cols[9].parse()?;
        if !s.gene_max_exons.contains_key(&gene_id) {
            s.gene_order.push(gene_id.clone());
            s.gene_max_exons.insert(gene_id.clone(), 0);
            s.gene_reads.insert(gene_id.clone(), Default::default());
            s.gene_trans.insert(gene_id.clone(), Default::default());
        }
        let e = s.gene_max_exons.get_mut(&gene_id).unwrap();
        if *e < num_exons {
            *e = num_exons;
        }
        s.gene_reads.get_mut(&gene_id).unwrap().insert(read_id);
        s.gene_trans.get_mut(&gene_id).unwrap().insert(trans_id.clone());
        let te = s.trans_max_exons.entry(trans_id).or_insert(0);
        if *te < num_exons {
            *te = num_exons;
        }
    }
    Ok(s)
}

/// Degradation signature statistics. Ports `tama_degradation_signature`.
fn degradation(
    capped: &std::path::Path,
    nocap: &std::path::Path,
    output: &std::path::Path,
) -> anyhow::Result<()> {
    let cap = read_deg(capped)?;
    let nc = read_deg(nocap)?;

    // (multi_exon_multi_read_trans_count, single_exon_gene, multi_exon_gene,
    //  single_exon_single_read_gene, multi_exon_single_read_gene)
    let tally = |s: &DegStats| -> (usize, usize, usize, usize, usize) {
        let (mut trans_count, mut se_gene, mut me_gene, mut se_sr_gene, mut me_sr_gene) =
            (0, 0, 0, 0, 0);
        for g in &s.gene_order {
            let reads = s.gene_reads[g].len();
            if s.gene_max_exons[g] < 2 {
                se_gene += 1;
                if reads == 1 {
                    se_sr_gene += 1;
                }
                continue;
            }
            me_gene += 1;
            if reads < 2 {
                me_sr_gene += 1;
                continue;
            }
            trans_count += s.gene_trans[g].len();
        }
        (trans_count, se_gene, me_gene, se_sr_gene, me_sr_gene)
    };
    let (cap_tc, cap_se_g, cap_me_g, cap_se_sr, cap_me_sr) = tally(&cap);
    let (nc_tc, nc_se_g, nc_me_g, nc_se_sr, nc_me_sr) = tally(&nc);

    let deg_sig = (cap_tc as f64 - nc_tc as f64) / cap_tc as f64;

    let se_trans = |s: &DegStats| s.trans_max_exons.values().filter(|&&e| e == 1).count();
    let me_trans = |s: &DegStats| s.trans_max_exons.values().filter(|&&e| e != 1).count();

    let mut out = tama_io::create_writer(output)?;
    let mut w = |line: String| writeln!(out, "{line}");
    w(format!("Degradation Signature = {}", fmt_float(deg_sig)))?;
    w(format!("Capped multi-exon, multi-read, transcript count = {cap_tc}"))?;
    w(format!("No-cap multi-exon, multi-read, transcript count = {nc_tc}"))?;
    w(format!("Capped total transcript count = {}", cap.trans_max_exons.len()))?;
    w(format!("No-cap total transcript count = {}", nc.trans_max_exons.len()))?;
    w(format!("Capped single exon trans count = {}", se_trans(&cap)))?;
    w(format!("No-cap single exon trans count = {}", se_trans(&nc)))?;
    w(format!("Capped multi exon trans count = {}", me_trans(&cap)))?;
    w(format!("No-cap multi exon trans count = {}", me_trans(&nc)))?;
    w(format!("Capped total gene count = {}", cap.gene_order.len()))?;
    w(format!("No-cap total gene count = {}", nc.gene_order.len()))?;
    w(format!("Capped single exon gene count = {cap_se_g}"))?;
    w(format!("No-cap single exon gene count = {nc_se_g}"))?;
    w(format!("Capped multi exon gene count = {cap_me_g}"))?;
    w(format!("No-cap multi exon gene count = {nc_me_g}"))?;
    w(format!("Capped single exon single read gene count = {cap_se_sr}"))?;
    w(format!("No-cap single exon single read gene count = {nc_se_sr}"))?;
    w(format!("Capped multi exon single read gene count = {cap_me_sr}"))?;
    w(format!("No-cap multi exon single read gene count = {nc_me_sr}"))?;
    Ok(())
}

/// Format a float like Python's `str()` (shortest round-trip; integers keep `.0`).
fn fmt_float(x: f64) -> String {
    if x == x.trunc() && x.is_finite() {
        format!("{x:.1}")
    } else {
        format!("{x}")
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
        if bin != 0 && read_count.is_multiple_of(bin) {
            writeln!(out, "{read_count}\t{gene_count}")?;
        }
        let _ = read_id;
        if seen.insert(gene.as_str()) {
            gene_count += 1;
        }
    }
    Ok(())
}
