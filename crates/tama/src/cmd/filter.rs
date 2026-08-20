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
    Fragments {
        /// BED file. (`-f`)
        #[arg(short = 'f', long)]
        bed: std::path::PathBuf,
        /// Output prefix. (`-o`)
        #[arg(short = 'o', long)]
        output: String,
        /// Exon/splice-junction wobble threshold. (`-m`)
        #[arg(short = 'm', long, default_value_t = 10)]
        wobble: i64,
        /// Transcript ends wobble threshold. (`-e`)
        #[arg(short = 'e', long, default_value_t = 500)]
        ends_wobble: i64,
        /// Single-exon overlap percent threshold. (`-s`)
        #[arg(short = 's', long, default_value_t = 20)]
        overlap_percent: i64,
    },
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
        Cmd::Fragments { bed, output, wobble, ends_wobble, overlap_percent } => {
            fragments(&bed, &output, wobble, ends_wobble, overlap_percent)
        }
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

/// A transcript for fragment filtering (mutable exon bounds during absorption).
#[derive(Clone)]
struct FragTx {
    scaffold: String,
    trans_start: i64,
    trans_end: i64,
    gene_id: String,
    trans_id: String,
    strand: char,
    num_exons: usize,
    exon_start_list: Vec<i64>,
    exon_end_list: Vec<i64>,
}

impl FragTx {
    fn start(&self) -> i64 {
        self.exon_start_list[0]
    }
    fn end(&self) -> i64 {
        *self.exon_end_list.last().unwrap()
    }
    fn length(&self) -> i64 {
        self.end() - self.start()
    }
    fn format_bed_line(&self) -> String {
        let t_start = self.trans_start.min(self.start());
        let t_end = self.trans_end.max(self.end());
        let (mut blocks, mut rel) = (String::new(), String::new());
        for k in 0..self.num_exons {
            if k > 0 {
                blocks.push(',');
                rel.push(',');
            }
            blocks.push_str(&(self.exon_end_list[k] - self.exon_start_list[k]).to_string());
            rel.push_str(&(self.exon_start_list[k] - t_start).to_string());
        }
        format!(
            "{}\t{t_start}\t{t_end}\t{};{}\t40\t{}\t{t_start}\t{t_end}\t255,0,0\t{}\t{blocks}\t{rel}",
            self.scaffold, self.gene_id, self.trans_id, self.strand, self.num_exons
        )
    }
}

/// Compare two transcripts for fragment absorption. Faithful port of
/// `compare_absorb_transcripts` (tama_cds/tama_id defaults). Returns
/// `Some((long_is_a, new_start0, new_end_last))` when `b`/`a` is a fragment of
/// the other, with the extended bounds to apply to the long model.
fn compare_absorb(
    a: &FragTx,
    b: &FragTx,
    wob: i64,
    ends_wob: i64,
    ovl_pct: i64,
) -> Option<(bool, i64, i64)> {
    let (a_ne, b_ne) = (a.num_exons, b.num_exons);
    let long_is_a = if a_ne == b_ne {
        if a.length() > b.length() {
            true
        } else if a.length() < b.length() {
            false
        } else {
            a.start() <= b.start()
        }
    } else {
        a_ne > b_ne
    };
    let (long, short) = if long_is_a { (a, b) } else { (b, a) };
    let (long_ne, short_ne) = (long.num_exons, short.num_exons);

    let mut new_start0 = long.start();
    let mut new_end_last = long.end();
    let mut matched;

    if long_ne == 1 && short_ne == 1 {
        let overlap = if long.start() <= short.start() {
            long.end() - short.start()
        } else {
            short.end() - long.start()
        };
        if overlap < 0 {
            return None;
        }
        let long_pct = overlap * 100 / long.length();
        let short_pct = overlap * 100 / short.length();
        if long_pct > ovl_pct || short_pct > ovl_pct {
            matched = true;
            new_start0 = long.start().min(short.start());
            new_end_last = long.end().max(short.end());
        } else {
            return None;
        }
    } else if short_ne == 1 {
        matched = false;
        for ei in 0..long_ne {
            let this_start = long.exon_start_list[ei];
            let this_end = long.exon_end_list[ei];
            if short.start() < this_end && short.end() > this_start {
                let start_wobble = this_start - short.start();
                let end_wobble = short.end() - this_end;
                if ei == 0 {
                    if start_wobble < ends_wob && end_wobble < wob {
                        matched = true;
                        if start_wobble > 0 {
                            new_start0 = short.start();
                        }
                    }
                } else if ei == long_ne - 1 {
                    if start_wobble < wob && end_wobble < ends_wob {
                        matched = true;
                        if end_wobble > 0 {
                            new_end_last = short.end();
                        }
                    }
                } else if start_wobble < wob && end_wobble < wob {
                    matched = true;
                }
            }
        }
        if !matched {
            return None;
        }
    } else {
        // both multi-exon: pairwise splice-junction matching
        let mut ls: std::collections::BTreeSet<usize> = Default::default();
        let mut le: std::collections::BTreeSet<usize> = Default::default();
        let mut ss: std::collections::BTreeSet<usize> = Default::default();
        let mut se: std::collections::BTreeSet<usize> = Default::default();
        for i in 0..long_ne {
            for j in 0..short_ne {
                let (lst, let_) = (long.exon_start_list[i], long.exon_end_list[i]);
                let (sst, set_) = (short.exon_start_list[j], short.exon_end_list[j]);
                if lst <= set_ && let_ >= sst {
                    let (sw, ew, sw_th, ew_th, ew_abs);
                    if i == 0 && j == 0 {
                        sw = lst - sst;
                        sw_th = ends_wob;
                        ew_abs = true;
                        ew = (set_ - let_).abs();
                        ew_th = wob;
                    } else if i > 0 && j == 0 {
                        sw = lst - sst;
                        sw_th = wob;
                        ew_abs = true;
                        ew = (set_ - let_).abs();
                        ew_th = wob;
                    } else if i == long_ne - 1 && j == short_ne - 1 {
                        sw = (lst - sst).abs();
                        sw_th = wob;
                        ew_abs = false;
                        ew = set_ - let_;
                        ew_th = ends_wob;
                    } else if i < long_ne - 1 && j == short_ne - 1 {
                        sw = (lst - sst).abs();
                        sw_th = wob;
                        ew_abs = false;
                        ew = set_ - let_;
                        ew_th = wob;
                    } else {
                        sw = (lst - sst).abs();
                        sw_th = wob;
                        ew_abs = true;
                        ew = (set_ - let_).abs();
                        ew_th = wob;
                    }
                    let _ = ew_abs;
                    if sw <= sw_th {
                        ls.insert(i);
                        ss.insert(j);
                    }
                    if ew <= ew_th {
                        le.insert(i);
                        se.insert(j);
                    }
                }
            }
        }
        let ls: Vec<usize> = ls.into_iter().collect();
        let le: Vec<usize> = le.into_iter().collect();
        let ss: Vec<usize> = ss.into_iter().collect();
        let se: Vec<usize> = se.into_iter().collect();
        let mut no_match = false;
        if ls.len() != le.len() || ss.len() != se.len() || ls.is_empty() || ss.is_empty() {
            no_match = true;
        }
        if ss.len() != short_ne || se.len() != short_ne {
            no_match = true;
        }
        if !no_match {
            // consecutive-index consistency (start==end index, +1 steps)
            for (idx, (&s, &e)) in ls.iter().zip(&le).enumerate() {
                if s != e {
                    no_match = true;
                }
                if idx > 0 && (s != ls[idx - 1] + 1 || e != le[idx - 1] + 1) {
                    no_match = true;
                }
            }
            for (idx, (&s, &e)) in ss.iter().zip(&se).enumerate() {
                if s != e {
                    no_match = true;
                }
                if idx > 0 && (s != ss[idx - 1] + 1 || e != se[idx - 1] + 1) {
                    no_match = true;
                }
            }
        }
        if !no_match && long_ne > short_ne {
            if ss.len() < short_ne {
                no_match = true;
            }
        } else if !no_match && (ls.len() < long_ne || ss.len() < short_ne) {
            no_match = true;
        }
        if no_match {
            return None;
        }
        matched = true;
        new_start0 = long.start().min(short.start());
        new_end_last = long.end().max(short.end());
    }

    if matched {
        Some((long_is_a, new_start0, new_end_last))
    } else {
        None
    }
}

/// Remove fragment models. Ports `tama_remove_fragment_models` (tama_id/tama_cds).
fn fragments(
    bed: &std::path::Path,
    output_prefix: &str,
    wobble: i64,
    ends_wobble: i64,
    overlap_percent: i64,
) -> anyhow::Result<()> {
    let mut gene_list: Vec<String> = Vec::new();
    let mut gene_trans: IndexMap<String, Vec<FragTx>> = IndexMap::new();
    for line in read_lines(bed)? {
        let cols: Vec<&str> = line.split('\t').collect();
        let bt = tama_core::model::BedTranscript::parse_bed_line(&line)?;
        let id_split: Vec<&str> = cols[3].split(';').collect();
        let gene_id = id_split[0].to_string();
        let trans_id = id_split[1].to_string();
        let tx = FragTx {
            scaffold: bt.scaffold.clone(),
            trans_start: bt.trans_start,
            trans_end: bt.trans_end,
            gene_id: gene_id.clone(),
            trans_id,
            strand: bt.strand,
            num_exons: bt.num_exons(),
            exon_start_list: bt.exon_start_list,
            exon_end_list: bt.exon_end_list,
        };
        if !gene_trans.contains_key(&gene_id) {
            gene_list.push(gene_id.clone());
        }
        gene_trans.entry(gene_id).or_default().push(tx);
    }

    let mut out_bed = tama_io::create_writer(format!("{output_prefix}.bed"))?;
    let mut out_discard = tama_io::create_writer(format!("{output_prefix}_discarded.txt"))?;

    for gene_id in &gene_list {
        let txs = gene_trans.get_mut(gene_id).unwrap();
        let n = txs.len();
        let id_to_idx: IndexMap<String, usize> =
            txs.iter().enumerate().map(|(i, t)| (t.trans_id.clone(), i)).collect();
        let mut removed: std::collections::HashSet<usize> = Default::default();

        for ai in 0..n {
            for bi in 0..n {
                if ai == bi || removed.contains(&ai) || removed.contains(&bi) {
                    continue;
                }
                if let Some((long_is_a, new_s0, new_el)) =
                    compare_absorb(&txs[ai], &txs[bi], wobble, ends_wobble, overlap_percent)
                {
                    let (long_idx, short_idx) = if long_is_a { (ai, bi) } else { (bi, ai) };
                    removed.insert(short_idx);
                    let long = &mut txs[long_idx];
                    let last = long.exon_end_list.len() - 1;
                    long.exon_start_list[0] = new_s0;
                    long.exon_end_list[last] = new_el;
                }
            }
        }
        let _ = &id_to_idx;

        for (i, tx) in txs.iter().enumerate() {
            if removed.contains(&i) {
                writeln!(out_discard, "{}", tx.format_bed_line())?;
            } else {
                writeln!(out_bed, "{}", tx.format_bed_line())?;
            }
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
