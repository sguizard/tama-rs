//! `tama filter` — transcript model filters.

use std::io::Write;

use anyhow::bail;
use clap::{Args as ClapArgs, Subcommand};
use indexmap::{IndexMap, IndexSet};

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
    SingleRead {
        /// Annotation BED file. (`-b`)
        #[arg(short = 'b', long)]
        bed: std::path::PathBuf,
        /// Read support levels file. (`-r`)
        #[arg(short = 'r', long)]
        read: std::path::PathBuf,
        /// Output prefix. (`-o`)
        #[arg(short = 'o', long)]
        output: String,
        /// Removal level: `gene` or `transcript`. (`-l`)
        #[arg(short = 'l', long, default_value = "gene")]
        level: String,
        /// Multi-exon handling: `keep_multi` or `remove_multi`. (`-k`)
        #[arg(short = 'k', long, default_value = "keep_multi")]
        multi: String,
        /// Minimum number of supporting sources. (`-s`)
        #[arg(short = 's', long, default_value_t = 1)]
        source_support: usize,
        /// Minimum number of supporting reads. (`-n`)
        #[arg(short = 'n', long, default_value_t = 2)]
        read_support: usize,
    },
}

pub fn run(args: Args) -> anyhow::Result<()> {
    match args.cmd {
        Cmd::SingleRead {
            bed, read, output, level, multi, source_support, read_support,
        } => single_read(&bed, &read, &output, &level, &multi, source_support, read_support),
        Cmd::PrimaryOrf => Err(super::not_implemented("filter primary-orf")),
        Cmd::Fragments => Err(super::not_implemented("filter fragments")),
        Cmd::Polya => Err(super::not_implemented("filter polya")),
    }
}

/// Remove single-read transcript models. Ports `tama_remove_single_read_models_levels`.
#[allow(clippy::too_many_arguments)]
fn single_read(
    bed: &std::path::Path,
    read: &std::path::Path,
    output_prefix: &str,
    level: &str,
    multi: &str,
    source_support: usize,
    read_support: usize,
) -> anyhow::Result<()> {
    // parse read support levels: trans -> reads; trans -> source -> reads
    let mut trans_reads: IndexMap<String, IndexSet<String>> = IndexMap::new();
    let mut source_reads: IndexMap<String, IndexMap<String, IndexSet<String>>> = IndexMap::new();
    for line in read_lines(read)? {
        if line.starts_with("merge_gene_id") {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        let trans_id = cols[1].to_string();
        trans_reads.entry(trans_id.clone()).or_default();
        let smap = source_reads.entry(trans_id.clone()).or_default();
        for src_read in cols[5].split(';') {
            let (src, reads) = src_read.split_once(':').unwrap_or((src_read, ""));
            let slot = smap.entry(src.to_string()).or_default();
            for r in reads.split(',') {
                trans_reads.get_mut(&trans_id).unwrap().insert(r.to_string());
                slot.insert(r.to_string());
            }
        }
    }

    // parse bed: gene order, per-gene trans, cols, num_exons
    let mut gene_list: Vec<String> = Vec::new();
    let mut gene_trans_list: IndexMap<String, Vec<String>> = IndexMap::new();
    let mut trans_cols: IndexMap<String, Vec<String>> = IndexMap::new();
    let mut trans_exons: IndexMap<String, usize> = IndexMap::new();
    for line in read_lines(bed)? {
        let cols: Vec<String> = line.split('\t').map(String::from).collect();
        let id_split: Vec<&str> = cols[3].split(';').collect();
        let (gene_id, trans_id) = (id_split[0].to_string(), id_split[1].to_string());
        let num_exons: usize = cols[9].parse()?;
        if !gene_trans_list.contains_key(&gene_id) {
            gene_list.push(gene_id.clone());
        }
        gene_trans_list.entry(gene_id).or_default().push(trans_id.clone());
        trans_exons.insert(trans_id.clone(), num_exons);
        trans_cols.insert(trans_id, cols);
    }

    let mut out_bed = tama_io::create_writer(format!("{output_prefix}.bed"))?;
    let mut out_report = tama_io::create_writer(format!("{output_prefix}_singleton_report.txt"))?;
    let mut out_single = tama_io::create_writer(format!("{output_prefix}_singleton.bed"))?;
    writeln!(out_report, "old_gene_id\told_trans_id\tsource_line\tnum_reads\tnew_gene_id\tnew_trans_id\tnum_exons")?;

    let source_line = |t: &str| -> String {
        source_reads.get(t).map(|m| m.keys().cloned().collect::<Vec<_>>().join(",")).unwrap_or_default()
    };
    let total_reads = |t: &str| -> usize { trans_reads.get(t).map(|s| s.len()).unwrap_or(0) };
    let num_sources = |t: &str| -> usize { source_reads.get(t).map(|m| m.len()).unwrap_or(0) };

    let mut new_gene_num = 0usize;

    for gene_id in &gene_list {
        let trans_list = gene_trans_list[gene_id].clone();

        // gene-level single-read removal for solo-transcript genes
        if trans_list.len() == 1 {
            let t = &trans_list[0];
            let num_exons = trans_exons[t];
            if num_sources(t) == 1 {
                let src = source_reads[t].keys().next().unwrap().clone();
                let num_reads = source_reads[t][&src].len();
                if num_reads == 1 && !(multi == "keep_multi" && num_exons > 1) {
                    writeln!(out_report, "{gene_id}\t{t}\t{}\t{num_reads}\tremoved_gene\tremoved_transcript\t{num_exons}", source_line(t))?;
                    writeln!(out_single, "{}", trans_cols[t].join("\t"))?;
                    continue;
                }
            }
        }

        let mut new_trans_num = 0usize;

        if level == "transcript" {
            let mut keep: IndexSet<String> = IndexSet::new();
            for t in &trans_list {
                let num_exons = trans_exons[t];
                let ns = num_sources(t);
                let tr = total_reads(t);
                // source_support==1 case: the Python's `ns>1 && tr>=thr` and
                // `tr>=thr` branches both reduce to `tr>=thr`.
                let keep_flag = if source_support > 1 {
                    ns >= source_support || (multi == "keep_multi" && num_exons > 1)
                } else {
                    tr >= read_support || (multi == "keep_multi" && num_exons > 1)
                };
                if keep_flag {
                    keep.insert(t.clone());
                }
            }
            let new_gene_id = if !keep.is_empty() {
                new_gene_num += 1;
                format!("G{new_gene_num}")
            } else {
                "removed_trans_level".to_string()
            };
            let mut report_lines: Vec<(String, bool)> = Vec::new();
            let all_removed = keep.is_empty();
            for t in &trans_list {
                let num_exons = trans_exons[t];
                let tr = total_reads(t);
                if keep.contains(t) {
                    new_trans_num += 1;
                    let new_trans_id = format!("{new_gene_id}.{new_trans_num}");
                    let mut cols = trans_cols[t].clone();
                    cols[3] = format!("{new_gene_id};{new_trans_id}");
                    writeln!(out_bed, "{}", cols.join("\t"))?;
                    report_lines.push((format!("{gene_id}\t{t}\t{}\t{tr}\t{new_gene_id}\t{new_trans_id}\t{num_exons}", source_line(t)), false));
                } else {
                    writeln!(out_single, "{}", trans_cols[t].join("\t"))?;
                    report_lines.push((format!("{gene_id}\t{t}\t{}\t{tr}\t{new_gene_id}\tremoved_transcript\t{num_exons}", source_line(t)), true));
                }
            }
            for (rl, _) in &report_lines {
                if all_removed {
                    // replace new_gene_id column (index 4) with all_trans_removed
                    let mut parts: Vec<&str> = rl.split('\t').collect();
                    parts[4] = "all_trans_removed";
                    writeln!(out_report, "{}", parts.join("\t"))?;
                } else {
                    writeln!(out_report, "{rl}")?;
                }
            }
        } else if level == "gene" {
            for t in &trans_list {
                let num_exons = trans_exons[t];
                let tr = total_reads(t);
                new_trans_num += 1;
                if new_trans_num == 1 {
                    new_gene_num += 1;
                }
                let new_gene_id = format!("G{new_gene_num}");
                let new_trans_id = format!("{new_gene_id}.{new_trans_num}");
                let mut cols = trans_cols[t].clone();
                cols[3] = format!("{new_gene_id};{new_trans_id}");
                writeln!(out_bed, "{}", cols.join("\t"))?;
                writeln!(out_report, "{gene_id}\t{t}\t{}\t{tr}\t{new_gene_id}\t{new_trans_id}\t{num_exons}", source_line(t))?;
            }
        } else {
            bail!("invalid -l level {level:?}");
        }
    }
    Ok(())
}

fn read_lines(path: &std::path::Path) -> anyhow::Result<Vec<String>> {
    let reader = tama_io::open_reader(path)?;
    let mut out = Vec::new();
    for line in std::io::BufRead::lines(reader) {
        let line = line?;
        if !line.trim().is_empty() {
            out.push(line);
        }
    }
    Ok(out)
}
