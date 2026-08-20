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
    /// Find model changes between two sources. (tama_find_model_changes)
    ModelChanges {
        /// Annotation BED file. (`-b`)
        #[arg(short = 'b', long)]
        bed: std::path::PathBuf,
        /// Read support levels file. (`-r`)
        #[arg(short = 'r', long)]
        read: std::path::PathBuf,
        /// Output prefix. (`-o`)
        #[arg(short = 'o', long)]
        output: String,
        /// Reference source name. (`-ref`)
        #[arg(long = "ref", default_value = "NA")]
        ref_source: String,
        /// Alternative source name. (`-alt`)
        #[arg(long = "alt", default_value = "NA")]
        alt_source: String,
    },
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
        Cmd::Degradation {
            capped,
            nocap,
            output,
        } => degradation(&capped, &nocap, &output),
        Cmd::ModelChanges {
            bed,
            read,
            output,
            ref_source,
            alt_source,
        } => model_changes(&bed, &read, &output, &ref_source, &alt_source),
        Cmd::Saturation {
            report,
            bin,
            output,
        } => saturation(&report, bin, &output),
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
        s.gene_trans
            .get_mut(&gene_id)
            .unwrap()
            .insert(trans_id.clone());
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
    w(format!(
        "Capped multi-exon, multi-read, transcript count = {cap_tc}"
    ))?;
    w(format!(
        "No-cap multi-exon, multi-read, transcript count = {nc_tc}"
    ))?;
    w(format!(
        "Capped total transcript count = {}",
        cap.trans_max_exons.len()
    ))?;
    w(format!(
        "No-cap total transcript count = {}",
        nc.trans_max_exons.len()
    ))?;
    w(format!(
        "Capped single exon trans count = {}",
        se_trans(&cap)
    ))?;
    w(format!(
        "No-cap single exon trans count = {}",
        se_trans(&nc)
    ))?;
    w(format!(
        "Capped multi exon trans count = {}",
        me_trans(&cap)
    ))?;
    w(format!("No-cap multi exon trans count = {}", me_trans(&nc)))?;
    w(format!(
        "Capped total gene count = {}",
        cap.gene_order.len()
    ))?;
    w(format!("No-cap total gene count = {}", nc.gene_order.len()))?;
    w(format!("Capped single exon gene count = {cap_se_g}"))?;
    w(format!("No-cap single exon gene count = {nc_se_g}"))?;
    w(format!("Capped multi exon gene count = {cap_me_g}"))?;
    w(format!("No-cap multi exon gene count = {nc_me_g}"))?;
    w(format!(
        "Capped single exon single read gene count = {cap_se_sr}"
    ))?;
    w(format!(
        "No-cap single exon single read gene count = {nc_se_sr}"
    ))?;
    w(format!(
        "Capped multi exon single read gene count = {cap_me_sr}"
    ))?;
    w(format!(
        "No-cap multi exon single read gene count = {nc_me_sr}"
    ))?;
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

type StrSet = IndexMap<String, ()>;
type NestedSet = IndexMap<String, StrSet>;

/// Find reads mapping to different genes/transcripts between two sources.
/// Ports `tama_find_model_changes`.
fn model_changes(
    bed: &std::path::Path,
    read: &std::path::Path,
    output_prefix: &str,
    ref_source: &str,
    alt_source: &str,
) -> anyhow::Result<()> {
    use std::io::BufRead;
    let read_lines = |p: &std::path::Path| -> anyhow::Result<Vec<String>> {
        let mut v = Vec::new();
        for l in tama_io::open_reader(p)?.lines() {
            let l = l?;
            if !l.trim().is_empty() {
                v.push(l);
            }
        }
        Ok(v)
    };

    // gene -> [chrom, min_start, max_end]
    let mut gene_pos: IndexMap<String, (String, i64, i64)> = IndexMap::new();
    for line in read_lines(bed)? {
        let c: Vec<&str> = line.split('\t').collect();
        let gene_id = c[3].split(';').next().unwrap_or("").to_string();
        let (chrom, ts, te): (String, i64, i64) = (c[0].to_string(), c[1].parse()?, c[2].parse()?);
        gene_pos
            .entry(gene_id)
            .and_modify(|g| {
                g.1 = g.1.min(ts);
                g.2 = g.2.max(te);
            })
            .or_insert((chrom, ts, te));
    }

    // read -> source -> {gene}; read -> {gene}; read -> source -> {trans}
    let mut read_src_gene: IndexMap<String, NestedSet> = IndexMap::new();
    let mut read_gene: IndexMap<String, StrSet> = IndexMap::new();
    let mut read_src_trans: IndexMap<String, NestedSet> = IndexMap::new();
    let mut gene_source: IndexMap<String, StrSet> = IndexMap::new();
    let mut trans_source: IndexMap<String, StrSet> = IndexMap::new();
    let mut source_set: StrSet = IndexMap::new();
    let mut all_reads: Vec<String> = Vec::new();

    for line in read_lines(read)? {
        if line.starts_with("merge_gene_id") {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        let (mg, mt, sources, support) = (c[0], c[1], c[4], c[5]);
        gene_source.entry(mg.to_string()).or_default();
        trans_source.entry(mt.to_string()).or_default();
        for s in sources.split(',') {
            gene_source.get_mut(mg).unwrap().insert(s.to_string(), ());
            trans_source.get_mut(mt).unwrap().insert(s.to_string(), ());
            source_set.insert(s.to_string(), ());
        }
        for src_read in support.split(';') {
            let (src, reads) = src_read.split_once(':').unwrap_or((src_read, ""));
            for r in reads.split(',') {
                if !read_src_gene.contains_key(r) {
                    all_reads.push(r.to_string());
                }
                read_src_gene
                    .entry(r.to_string())
                    .or_default()
                    .entry(src.to_string())
                    .or_default()
                    .insert(mg.to_string(), ());
                read_gene
                    .entry(r.to_string())
                    .or_default()
                    .insert(mg.to_string(), ());
                read_src_trans
                    .entry(r.to_string())
                    .or_default()
                    .entry(src.to_string())
                    .or_default()
                    .insert(mt.to_string(), ());
            }
        }
    }

    let mut out_gene = tama_io::create_writer(format!("{output_prefix}_diff_genes.txt"))?;
    let mut out_trans = tama_io::create_writer(format!("{output_prefix}_diff_trans.txt"))?;
    let mut out_report = tama_io::create_writer(format!("{output_prefix}_diff_report.txt"))?;
    let mut out_onegene =
        tama_io::create_writer(format!("{output_prefix}_diff_one_source_genes.txt"))?;
    let mut out_onetrans =
        tama_io::create_writer(format!("{output_prefix}_diff_one_source_trans.txt"))?;
    writeln!(
        out_gene,
        "read_id\tnum_genes\tall_gene_line\tall_pos_line\tall_trans_line"
    )?;
    writeln!(out_trans, "read_id\talt_trans_diff_count\talt_diff_trans_id_list_line\talt_trans_id_list_line\tref_trans_id_list_line")?;

    let mut all_source_list: Vec<String> = source_set.keys().cloned().collect();
    all_source_list.sort();
    let mut diff_gene_src: IndexMap<String, StrSet> = IndexMap::new();
    let mut diff_trans_src: IndexMap<String, StrSet> = IndexMap::new();
    for s in &all_source_list {
        diff_gene_src.insert(s.clone(), IndexMap::new());
        diff_trans_src.insert(s.clone(), IndexMap::new());
    }
    let mut merge_diff_gene: StrSet = IndexMap::new();
    let mut merge_diff_trans: StrSet = IndexMap::new();
    let mut read_diff_gene: StrSet = IndexMap::new();
    let mut read_diff_trans_count = 0usize;

    for read_id in &all_reads {
        let num_genes = read_gene[read_id].len();
        if num_genes > 1 {
            let (mut gene_lines, mut pos_lines, mut trans_lines) =
                (Vec::new(), Vec::new(), Vec::new());
            for (src, genes) in &read_src_gene[read_id] {
                for g in genes.keys() {
                    let (chrom, gs, ge) = &gene_pos[g];
                    pos_lines.push(format!("{chrom}:{gs}-{ge}"));
                    gene_lines.push(format!("{src}:{g}"));
                    diff_gene_src.get_mut(src).unwrap().insert(g.clone(), ());
                    merge_diff_gene.insert(g.clone(), ());
                }
            }
            for (src, ts) in &read_src_trans[read_id] {
                for t in ts.keys() {
                    trans_lines.push(format!("{src}:{t}"));
                }
            }
            read_diff_gene.insert(read_id.clone(), ());
            writeln!(
                out_gene,
                "{read_id}\t{num_genes}\t{}\t{}\t{}",
                gene_lines.join(","),
                pos_lines.join(","),
                trans_lines.join(",")
            )?;
            continue;
        }

        // transcript differences
        let empty = IndexMap::new();
        let src_trans = read_src_trans.get(read_id).unwrap_or(&empty);
        let ref_trans: StrSet = src_trans.get(ref_source).cloned().unwrap_or_default();
        let (mut alt_list, mut alt_diff, mut diff_count) = (Vec::new(), Vec::new(), 0usize);
        let mut diff_flag = false;
        if let Some(alt_ts) = src_trans.get(alt_source) {
            for t in alt_ts.keys() {
                alt_list.push(t.clone());
                if !ref_trans.contains_key(t) && src_trans.contains_key(ref_source) {
                    diff_flag = true;
                    diff_count += 1;
                    alt_diff.push(t.clone());
                    diff_trans_src
                        .get_mut(alt_source)
                        .unwrap()
                        .insert(t.clone(), ());
                    merge_diff_trans.insert(t.clone(), ());
                }
            }
        }
        if !diff_flag {
            continue;
        }
        let mut ref_list = Vec::new();
        if let Some(ref_ts) = src_trans.get(ref_source) {
            for t in ref_ts.keys() {
                ref_list.push(t.clone());
                diff_trans_src
                    .get_mut(ref_source)
                    .unwrap()
                    .insert(t.clone(), ());
                merge_diff_trans.insert(t.clone(), ());
            }
        } else {
            continue; // read discarded in ref
        }
        read_diff_trans_count += 1;
        writeln!(
            out_trans,
            "{read_id}\t{diff_count}\t{}\t{}\t{}",
            alt_diff.join(","),
            alt_list.join(","),
            ref_list.join(",")
        )?;
    }

    writeln!(out_report, "num_diff_gene_reads: {}", read_diff_gene.len())?;
    writeln!(out_report, "num_diff_trans_reads: {read_diff_trans_count}")?;
    writeln!(out_report, "num_merge_diff_gene: {}", merge_diff_gene.len())?;
    writeln!(
        out_report,
        "num_merge_diff_trans: {}",
        merge_diff_trans.len()
    )?;
    for s in &all_source_list {
        writeln!(
            out_report,
            "this_source_diff_genes {s}: {}",
            diff_gene_src[s].len()
        )?;
        writeln!(
            out_report,
            "this_source_diff_trans {s}: {}",
            diff_trans_src[s].len()
        )?;
    }

    writeln!(out_onegene, "merge_source\tmerge_gene_id")?;
    writeln!(out_onetrans, "merge_source\tmerge_trans_id")?;
    let mut only_gene: IndexMap<String, StrSet> = IndexMap::new();
    let mut only_trans: IndexMap<String, StrSet> = IndexMap::new();
    for s in &all_source_list {
        only_gene.insert(s.clone(), IndexMap::new());
        only_trans.insert(s.clone(), IndexMap::new());
    }
    for g in merge_diff_gene.keys() {
        let srcs: Vec<String> = gene_source[g].keys().cloned().collect();
        if srcs.len() == 1 {
            only_gene.get_mut(&srcs[0]).unwrap().insert(g.clone(), ());
            writeln!(out_onegene, "{}\t{g}", srcs[0])?;
        }
    }
    for t in merge_diff_trans.keys() {
        let srcs: Vec<String> = trans_source[t].keys().cloned().collect();
        if srcs.len() == 1 {
            let gene = t.split('.').next().unwrap_or(t);
            if !only_gene[&srcs[0]].contains_key(gene) {
                only_trans.get_mut(&srcs[0]).unwrap().insert(t.clone(), ());
                writeln!(out_onetrans, "{}\t{t}", srcs[0])?;
            }
        }
    }
    let (mut tot_g, mut tot_t) = (0usize, 0usize);
    for s in &all_source_list {
        tot_g += only_gene[s].len();
        tot_t += only_trans[s].len();
        writeln!(
            out_report,
            "only_source_num_genes {s}: {}",
            only_gene[s].len()
        )?;
        writeln!(
            out_report,
            "only_source_num_trans {s}: {}",
            only_trans[s].len()
        )?;
    }
    writeln!(out_report, "total_one_source_genes_count: {tot_g}")?;
    writeln!(out_report, "total_one_source_trans_count: {tot_t}")?;
    Ok(())
}

/// Sampling saturation curve. Ports `tama_sampling_saturation_curve`.
///
/// Note: the curve is cumulative over reads in the order they first appear; the
/// original iterated a Python-2 dict (hash order), so only the final totals are
/// guaranteed to match — the intermediate sampling points differ.
fn saturation(
    report: &std::path::Path,
    bin: usize,
    output: &std::path::Path,
) -> anyhow::Result<()> {
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
