//! `tama merge` — merge transcript annotations across sources.
//!
//! Ports `tama_merge.py`. Reads a filelist of BED sources (each with a cap flag
//! and start/junction/end priorities), groups overlapping transcripts into genes,
//! collapses matching models (priority-aware), and writes the merged BED plus
//! merge/trans/gene reports.
//!
//! Supports capped, no_cap, and mixed sources via the phased grouping in
//! `group_gene` (capped connected components, then no_cap attachment, then
//! no_cap-only components), using the capped / capped-nocap / both-nocap
//! comparisons from the original. The `-s` (source gene/trans ids) and `-cds`
//! (source CDS) options are supported.

use std::collections::BTreeMap;
use std::io::Write;

use anyhow::{bail, Context};
use clap::Parser;
use indexmap::{IndexMap, IndexSet};

use crate::cmd::opts::{Dup, EndsOpt};
use tama_core::gene::{gene_group, GeneMember};
use tama_core::model::BedTranscript;

#[derive(Parser)]
pub struct Args {
    /// File list describing the annotations to merge. (`-f`)
    #[arg(short = 'f', long = "filelist")]
    pub filelist: std::path::PathBuf,
    /// Output prefix. (`-p`)
    #[arg(short = 'p', long = "prefix")]
    pub prefix: String,
    /// Collapse exon ends. (`-e`)
    #[arg(short = 'e', long, value_enum, default_value = "common_ends")]
    pub ends: EndsOpt,
    /// 5' threshold. (`-a`, merge default 20)
    #[arg(short = 'a', long, default_value_t = 20)]
    pub five_prime: i64,
    /// Exon/splice-junction threshold. (`-m`)
    #[arg(short = 'm', long, default_value_t = 10)]
    pub exon_thresh: i64,
    /// 3' threshold. (`-z`, merge default 20)
    #[arg(short = 'z', long, default_value_t = 20)]
    pub three_prime: i64,
    /// Duplicate merge behaviour. (`-d`)
    #[arg(short = 'd', long, value_enum, default_value = "no_merge")]
    pub dup: Dup,
    /// Use gene/transcript IDs from this merge source. (`-s`)
    #[arg(short = 's', long)]
    pub source_id: Option<String>,
    /// Use CDS from this merge source. (`-cds`)
    #[arg(long = "cds")]
    pub cds_source: Option<String>,
}

/// source id -> (filename, seq_type, (start, junction, end) priorities)
type SourceMap = IndexMap<String, (String, String, (i64, i64, i64))>;

struct Thresholds {
    five: i64,
    three: i64,
    exon: i64,
    longest_ends: bool,
}

#[derive(Clone)]
struct MergeTx {
    uniq_trans_id: String,
    source_id: String,
    /// `true` for capped sources, `false` for no_cap.
    capped: bool,
    scaffold: String,
    strand: char,
    trans_start: i64,
    trans_end: i64,
    num_exons: usize,
    exon_start_list: Vec<i64>,
    exon_end_list: Vec<i64>,
    start_priority: i64,
    junction_priority: i64,
    end_priority: i64,
    /// Original (un-prefixed) gene/trans ids and CDS from the source BED — used by
    /// the `-s`/`-cds` options.
    src_gene_id: String,
    src_trans_id: String,
    src_cds_start: i64,
    src_cds_end: i64,
}

struct Collapsed {
    start: Vec<i64>,
    end: Vec<i64>,
    start_wobble: Vec<i64>,
    end_wobble: Vec<i64>,
    /// support: per collapsed exon, sorted uniq_trans_ids supporting the start/end.
    e_start_support: Vec<String>,
    e_end_support: Vec<String>,
}

struct MergedTrans {
    trans_id: String,
    scaffold: String,
    strand: char,
    num_exons: usize,
    collapsed: Collapsed,
    members: Vec<String>, // uniq_trans_ids in insertion order
    rgb: String,
}

impl MergedTrans {
    fn start_pos(&self) -> i64 {
        self.collapsed.start[0]
    }
    fn end_pos(&self) -> i64 {
        *self.collapsed.end.last().unwrap()
    }

    /// `extra_ids` are appended to the col-4 id line (`-s`); `cds` overrides the
    /// thick start/end (`-cds`).
    fn format_bed_line(&self, extra_ids: &[(String, String)], cds: Option<(i64, i64)>) -> String {
        let gene_id = self.trans_id.split('.').next().unwrap_or(&self.trans_id);
        let mut id_line = format!("{};{}", gene_id, self.trans_id);
        for (g, t) in extra_ids {
            id_line.push_str(&format!(";{g};{t}"));
        }
        let (thick_start, thick_end) = match cds {
            Some((s, e)) if s != 0 => (s, e),
            _ => (self.start_pos(), self.end_pos()),
        };
        let (mut sizes, mut starts) = (String::new(), String::new());
        for k in 0..self.num_exons {
            if k > 0 {
                sizes.push(',');
                starts.push(',');
            }
            sizes.push_str(&(self.collapsed.end[k] - self.collapsed.start[k]).to_string());
            starts.push_str(&(self.collapsed.start[k] - self.start_pos()).to_string());
        }
        [
            self.scaffold.clone(),
            self.start_pos().to_string(),
            self.end_pos().to_string(),
            id_line,
            "40".to_string(),
            self.strand.to_string(),
            thick_start.to_string(),
            thick_end.to_string(),
            self.rgb.clone(),
            self.num_exons.to_string(),
            sizes,
            starts,
        ]
        .join("\t")
    }

    /// Sort key: interleaved exon starts/ends only (no leading trans start/end),
    /// matching merge's `sort_transcripts`.
    fn pos_key(&self) -> Vec<i64> {
        let mut s = self.collapsed.start.clone();
        let mut e = self.collapsed.end.clone();
        s.sort_unstable();
        e.sort_unstable();
        let mut key = Vec::with_capacity(s.len() * 2);
        for k in 0..s.len() {
            key.push(s[k]);
            key.push(e[k]);
        }
        key
    }
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let th = Thresholds {
        five: args.five_prime,
        three: args.three_prime,
        exon: args.exon_thresh,
        longest_ends: args.ends == EndsOpt::LongestEnds,
    };
    // ---- parse filelist ----
    let filelist_dir = args
        .filelist
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_default();
    let filelist = std::fs::read_to_string(&args.filelist)
        .with_context(|| format!("reading filelist {}", args.filelist.display()))?;
    // source -> (filename, seq_type, priority_rank)
    let mut sources: SourceMap = IndexMap::new();
    for line in filelist.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 4 {
            bail!("filelist line must have 4 tab-separated fields: {line:?}");
        }
        let (filename, seq_type, prio, source_id) = (f[0], f[1], f[2], f[3]);
        if seq_type != "capped" && seq_type != "no_cap" {
            bail!("seq type must be capped or no_cap: {seq_type:?}");
        }
        let pr: Vec<i64> = prio
            .split(',')
            .filter_map(|x| x.trim().parse().ok())
            .collect();
        if pr.len() != 3 {
            bail!("priority rank must be start,junction,end: {prio:?}");
        }
        sources.entry(source_id.to_string()).or_insert((
            filename.to_string(),
            seq_type.to_string(),
            (pr[0], pr[1], pr[2]),
        ));
    }

    // ---- -s / -cds source selections (comma-separated source names) ----
    let parse_sources = |flag: &Option<String>, label: &str| -> anyhow::Result<IndexSet<String>> {
        let mut set = IndexSet::new();
        if let Some(v) = flag {
            for name in v.split(',').map(|s| s.trim()).filter(|s| !s.is_empty()) {
                if !sources.contains_key(name) {
                    bail!("{label} source {name:?} is not in the filelist");
                }
                set.insert(name.to_string());
            }
        }
        Ok(set)
    };
    let source_id_flags = parse_sources(&args.source_id, "-s")?;
    let cds_flags = parse_sources(&args.cds_source, "-cds")?;

    // source colour (sorted sources -> 1..=10)
    let mut colour_sources: Vec<&String> = sources.keys().collect();
    colour_sources.sort();
    let source_colour: IndexMap<String, usize> = colour_sources
        .iter()
        .enumerate()
        .map(|(i, s)| ((*s).clone(), (i + 1).min(10)))
        .collect();

    // ---- read beds into bed_dict[scaffold][start][end] = Vec<MergeTx> ----
    type EndMap = BTreeMap<i64, Vec<MergeTx>>;
    type StartMap = BTreeMap<i64, EndMap>;
    let mut bed_dict: BTreeMap<String, StartMap> = BTreeMap::new();

    let mut total_loaded = 0usize;
    for (source_id, (filename, seq, prio)) in &sources {
        let capped = seq == "capped";
        let path = if std::path::Path::new(filename).is_absolute() {
            std::path::PathBuf::from(filename)
        } else {
            filelist_dir.join(filename)
        };
        let reader = tama_io::open_reader(&path)
            .with_context(|| format!("opening bed {}", path.display()))?;
        let transcripts = tama_core::bed::read_bed(reader)
            .with_context(|| format!("parsing bed {}", path.display()))?;
        let n = transcripts.len();
        for bt in transcripts {
            let tx = to_merge_tx(&bt, source_id, capped, *prio);
            bed_dict
                .entry(tx.scaffold.clone())
                .or_default()
                .entry(tx.trans_start)
                .or_default()
                .entry(tx.trans_end)
                .or_default()
                .push(tx);
        }
        total_loaded += n;
        log::debug!("merge: loaded {n} transcripts from source {source_id}");
    }
    log::debug!(
        "merge: {total_loaded} transcripts from {} sources across {} scaffolds; grouping…",
        sources.len(),
        bed_dict.len()
    );

    // ---- output writers ----
    let p = &args.prefix;
    let mut out_bed = tama_io::create_writer(format!("{p}.bed"))?;
    let mut out_merge = tama_io::create_writer(format!("{p}_merge.txt"))?;
    let mut out_trans = tama_io::create_writer(format!("{p}_trans_report.txt"))?;
    let mut out_gene = tama_io::create_writer(format!("{p}_gene_report.txt"))?;
    writeln!(out_trans, "transcript_id\tnum_clusters\tsources\tstart_wobble_list\tend_wobble_list\texon_start_support\texon_end_support\tall_source_trans")?;
    writeln!(out_gene, "gene_id\tnum_clusters\tnum_final_trans\tsources\tchrom\tstart\tend\tsource_genes\tsource_summary")?;

    // ---- position grouping per scaffold, then process each group ----
    let mut total_gene_count = 0i64;
    let num_scaffolds = bed_dict.len();
    for (scaff_idx, (scaffold, start_map)) in bed_dict.iter().enumerate() {
        let starts: Vec<i64> = start_map.keys().copied().collect();
        let num = starts.len();
        let mut group: Vec<MergeTx> = Vec::new();
        let mut last_trans_end = 0i64;

        for (si, &trans_start) in starts.iter().enumerate() {
            let end_map = &start_map[&trans_start];
            let ends: Vec<i64> = end_map.keys().copied().collect();
            if last_trans_end == 0 {
                last_trans_end = ends[0];
            }
            if trans_start < last_trans_end {
                for &te in &ends {
                    group.extend(end_map[&te].iter().cloned());
                    if last_trans_end < te {
                        last_trans_end = te;
                    }
                }
            } else {
                total_gene_count = process_trans_group(
                    &group,
                    total_gene_count,
                    &source_colour,
                    &th,
                    &source_id_flags,
                    &cds_flags,
                    &mut out_bed,
                    &mut out_merge,
                    &mut out_trans,
                    &mut out_gene,
                )?;
                group = Vec::new();
                for &te in &ends {
                    group.extend(end_map[&te].iter().cloned());
                    if last_trans_end < te {
                        last_trans_end = te;
                    }
                }
            }
            if si + 1 == num {
                total_gene_count = process_trans_group(
                    &group,
                    total_gene_count,
                    &source_colour,
                    &th,
                    &source_id_flags,
                    &cds_flags,
                    &mut out_bed,
                    &mut out_merge,
                    &mut out_trans,
                    &mut out_gene,
                )?;
            }
        }
        log::debug!(
            "merge: scaffold {}/{} ({scaffold}) done — {total_gene_count} genes so far",
            scaff_idx + 1,
            num_scaffolds
        );
    }

    log::info!("merge done: {total_gene_count} genes");
    Ok(())
}

fn to_merge_tx(
    bt: &BedTranscript,
    source_id: &str,
    capped: bool,
    prio: (i64, i64, i64),
) -> MergeTx {
    MergeTx {
        uniq_trans_id: format!("{}_{}", source_id, bt.trans_id()),
        source_id: source_id.to_string(),
        capped,
        scaffold: bt.scaffold.clone(),
        strand: bt.strand,
        trans_start: bt.trans_start,
        trans_end: bt.trans_end,
        num_exons: bt.num_exons(),
        exon_start_list: bt.exon_start_list.clone(),
        exon_end_list: bt.exon_end_list.clone(),
        start_priority: prio.0,
        junction_priority: prio.1,
        end_priority: prio.2,
        src_gene_id: bt.gene_id().to_string(),
        src_trans_id: bt.trans_id().to_string(),
        src_cds_start: bt.cds_start,
        src_cds_end: bt.cds_end,
    }
}

#[allow(clippy::too_many_arguments)]
fn process_trans_group(
    group: &[MergeTx],
    mut total_gene_count: i64,
    source_colour: &IndexMap<String, usize>,
    th: &Thresholds,
    source_id_flags: &IndexSet<String>,
    cds_flags: &IndexSet<String>,
    out_bed: &mut Box<dyn Write>,
    out_merge: &mut Box<dyn Write>,
    out_trans: &mut Box<dyn Write>,
    out_gene: &mut Box<dyn Write>,
) -> anyhow::Result<i64> {
    if group.is_empty() {
        return Ok(total_gene_count);
    }
    // index by uniq_trans_id
    let by_id: IndexMap<String, MergeTx> = group
        .iter()
        .map(|t| (t.uniq_trans_id.clone(), t.clone()))
        .collect();

    // split by strand
    let fwd: Vec<&MergeTx> = group.iter().filter(|t| t.strand == '+').collect();
    let rev: Vec<&MergeTx> = group.iter().filter(|t| t.strand == '-').collect();

    let fwd_genes = gene_groups(&fwd);
    let rev_genes = gene_groups(&rev);
    let mut gene_entries: Vec<(i64, u8, Vec<String>)> = Vec::new();
    for g in fwd_genes {
        gene_entries.push((g.0, 0, g.1));
    }
    for g in rev_genes {
        gene_entries.push((g.0, 1, g.1));
    }
    gene_entries.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    for (_start, _rev, ids) in gene_entries {
        total_gene_count += 1;
        let gene_txs: Vec<MergeTx> = ids.iter().map(|id| by_id[id].clone()).collect();

        let match_groups = group_gene(&gene_txs, th);

        let mut merged: Vec<MergedTrans> = Vec::new();
        for idx_group in &match_groups {
            let members: Vec<MergeTx> = idx_group.iter().map(|&k| gene_txs[k].clone()).collect();
            let collapsed = collapse_transcripts(&members, th);
            let num_exons = members.iter().map(|t| t.num_exons).max().unwrap();
            let rgb = rgb_for(&members, source_colour);
            merged.push(MergedTrans {
                trans_id: String::new(),
                scaffold: members[0].scaffold.clone(),
                strand: members[0].strand,
                num_exons,
                collapsed,
                members: members.iter().map(|t| t.uniq_trans_id.clone()).collect(),
                rgb,
            });
        }

        sort_merged(&mut merged);

        let mut trans_count = 0i64;
        let mut track_gene_source: IndexSet<String> = IndexSet::new();
        for mut m in merged {
            trans_count += 1;
            m.trans_id = format!("G{total_gene_count}.{trans_count}");

            writeln!(out_trans, "{}", format_trans_report(&m))?;

            // -s: append the source's original gene;trans ids; -cds: override CDS.
            let mut extra_ids: Vec<(String, String)> = Vec::new();
            let mut cds: Option<(i64, i64)> = None;
            for uid in &m.members {
                // source gene id (uniq_trans_id minus the source prefix, gene part)
                let src_trans = uid.split_once('_').map(|x| x.1).unwrap_or(uid);
                let src_name = uid.split('_').next().unwrap_or("");
                let gene_part = src_trans.split('.').next().unwrap_or(src_trans);
                track_gene_source.insert(format!("{src_name}_{gene_part}"));

                let tx = &by_id[uid];
                if source_id_flags.contains(&tx.source_id) {
                    extra_ids.push((tx.src_gene_id.clone(), tx.src_trans_id.clone()));
                }
                if cds_flags.contains(&tx.source_id) {
                    cds = Some((tx.src_cds_start, tx.src_cds_end));
                }
                writeln!(out_merge, "{}", member_bed_line(tx, &m.trans_id))?;
            }

            writeln!(out_bed, "{}", m.format_bed_line(&extra_ids, cds))?;
        }

        let mut src_genes: Vec<String> = track_gene_source.into_iter().collect();
        src_genes.sort();
        writeln!(
            out_gene,
            "{}",
            format_gene_report(&gene_txs, total_gene_count, trans_count, &src_genes)
        )?;
    }

    Ok(total_gene_count)
}

/// Gene-group a strand's transcripts, returning (gene_start, member uniq ids).
fn gene_groups(txs: &[&MergeTx]) -> Vec<(i64, Vec<String>)> {
    if txs.is_empty() {
        return Vec::new();
    }
    let members: Vec<GeneMember> = txs
        .iter()
        .map(|t| GeneMember {
            id: &t.uniq_trans_id,
            exon_starts: &t.exon_start_list,
            exon_ends: &t.exon_end_list,
        })
        .collect();
    gene_group(&members)
        .into_iter()
        .map(|g| (g.gene_start, g.trans_ids))
        .collect()
}

/// Capped comparison: same exon count and every boundary matches within the
/// appropriate threshold. Mirrors `compare_transcripts_both_capped`.
fn same_transcript_capped(a: &MergeTx, b: &MergeTx, th: &Thresholds) -> bool {
    if a.num_exons != b.num_exons {
        return false;
    }
    let n = a.num_exons;
    let strand = a.strand;
    for i in 0..n {
        let (ja, jb) = if strand == '+' {
            (a.num_exons - 1 - i, b.num_exons - 1 - i)
        } else {
            (i, i)
        };
        let mut start_th = th.exon;
        let mut end_th = th.exon;
        if strand == '+' {
            if i == 0 {
                end_th = th.three;
            }
            if i == n - 1 {
                start_th = th.five;
            }
        } else {
            if i == 0 {
                start_th = th.three;
            }
            if i == n - 1 {
                end_th = th.five;
            }
        }
        if (a.exon_start_list[ja] - b.exon_start_list[jb]).abs() > start_th {
            return false;
        }
        if (a.exon_end_list[ja] - b.exon_end_list[jb]).abs() > end_th {
            return false;
        }
    }
    true
}

/// Order a pair into (long, short). Returns `None` when the pair can't match:
/// for a capped+no_cap pair the no_cap must not have more exons than the capped
/// (`compare_transcripts_capped_nocap`). For two no_cap transcripts the longer
/// (more exons; ties broken by the more-extended 5' end) is the "long" one.
fn assign_long_short<'a>(
    a: &'a MergeTx,
    b: &'a MergeTx,
    strand: char,
) -> Option<(&'a MergeTx, &'a MergeTx)> {
    if a.capped != b.capped {
        let (capped, nocap) = if a.capped { (a, b) } else { (b, a) };
        if nocap.num_exons > capped.num_exons {
            return None;
        }
        return Some((capped, nocap));
    }
    // both no_cap
    if a.num_exons != b.num_exons {
        if a.num_exons > b.num_exons {
            Some((a, b))
        } else {
            Some((b, a))
        }
    } else if strand == '+' {
        // long = the one whose 5' (first) exon start is earlier/equal
        if a.exon_start_list[0] <= b.exon_start_list[0] {
            Some((a, b))
        } else {
            Some((b, a))
        }
    } else if a.exon_end_list[a.num_exons - 1] >= b.exon_end_list[b.num_exons - 1] {
        Some((a, b))
    } else {
        Some((b, a))
    }
}

/// Comparison when at least one transcript is no_cap. Returns `true` when the
/// original would yield `same_transcript` / `same_three_prime_same_exons` /
/// `same_three_prime_diff_exons` (i.e. a 5'-degraded model joins its longer
/// relative). Faithful boolean reduction of `compare_transcripts_capped_nocap`
/// and `compare_transcripts_both_nocap`.
fn nocap_match(a: &MergeTx, b: &MergeTx, th: &Thresholds) -> bool {
    let strand = a.strand;
    let Some((long, short)) = assign_long_short(a, b, strand) else {
        return false;
    };
    let max = long.num_exons;
    let min = short.num_exons;
    if min == 0 {
        return false;
    }
    // Iterate from the 3' end (strand-corrected), comparing the shortest overlap.
    for i in 0..min {
        let (jl, js) = if strand == '+' {
            (max - 1 - i, min - 1 - i)
        } else {
            (i, i)
        };
        let mut start_th = th.exon;
        let mut end_th = th.exon;
        if strand == '+' {
            if i == 0 {
                end_th = th.three;
            }
            if i == max - 1 {
                start_th = th.five;
            }
        } else {
            if i == 0 {
                start_th = th.three;
            }
            if i == max - 1 {
                end_th = th.five;
            }
        }
        let ls = long.exon_start_list[jl];
        let ss = short.exon_start_list[js];
        let le = long.exon_end_list[jl];
        let se = short.exon_end_list[js];
        let start_match = (ls - ss).abs() <= start_th;
        let end_match = (le - se).abs() <= end_th;

        if i < min - 1 {
            // internal exon: both boundaries must match
            if !start_match || !end_match {
                return false;
            }
        } else {
            // 5'-most compared exon (i == min - 1): allow 5' degradation
            if strand == '+' {
                if !end_match {
                    return false;
                }
                if !start_match && ls > ss {
                    // long's 5' start is later than short's → not a degraded form
                    return false;
                }
            } else {
                if !start_match {
                    return false;
                }
                if !end_match && le < se {
                    return false;
                }
            }
        }
    }
    true
}

/// Group a gene's transcripts into merged models, mirroring the phased
/// `simplify_gene` in the original:
///
/// 1. **Capped groups** — connected components over `same_transcript_capped`
///    (`hunter_prey_capped`). Two distinct capped groups never merge.
/// 2. **Attach no_cap** — each no_cap transcript is added to *every* capped group
///    that contains a capped member it matches (`hunter_prey_mixed`); a no_cap can
///    thus support more than one merged model, and never bridges two capped groups.
/// 3. **No_cap-only groups** — no_caps not attached to any capped group are grouped
///    among themselves by connected components over `nocap_match`
///    (`hunter_prey_nocap`).
///
/// Returns member index lists (into `txs`); a no_cap index may appear in several.
fn group_gene(txs: &[MergeTx], th: &Thresholds) -> Vec<Vec<usize>> {
    let capped_idx: Vec<usize> = (0..txs.len()).filter(|&i| txs[i].capped).collect();
    let nocap_idx: Vec<usize> = (0..txs.len()).filter(|&i| !txs[i].capped).collect();

    // Phase 1: capped connected components.
    let mut groups = connected_components(&capped_idx, |a, b| {
        same_transcript_capped(&txs[a], &txs[b], th)
    });

    // Phase 2: attach no_caps to matching capped groups.
    let mut attached = vec![false; txs.len()];
    for g in groups.iter_mut() {
        for &n in &nocap_idx {
            let matches = g.iter().any(|&h| {
                txs[h].num_exons >= txs[n].num_exons && nocap_match(&txs[h], &txs[n], th)
            });
            if matches {
                g.push(n);
                attached[n] = true;
            }
        }
    }

    // Phase 3: remaining no_caps grouped among themselves.
    let remaining: Vec<usize> = nocap_idx.into_iter().filter(|&n| !attached[n]).collect();
    groups.extend(connected_components(&remaining, |a, b| {
        nocap_match(&txs[a], &txs[b], th)
    }));

    groups
}

/// Connected components over `idxs` under the symmetric `matches` relation.
/// Returns groups of original indices, first-seen order preserved.
fn connected_components(idxs: &[usize], matches: impl Fn(usize, usize) -> bool) -> Vec<Vec<usize>> {
    let n = idxs.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        let mut c = x;
        while p[c] != r {
            let nx = p[c];
            p[c] = r;
            c = nx;
        }
        r
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if matches(idxs[i], idxs[j]) {
                let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }
    let mut groups: IndexMap<usize, Vec<usize>> = IndexMap::new();
    for (i, &orig) in idxs.iter().enumerate() {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(orig);
    }
    groups.into_values().collect()
}

/// Priority-aware exon collapsing. Faithful port of the merge `collapse_transcripts`.
fn collapse_transcripts(txs: &[MergeTx], th: &Thresholds) -> Collapsed {
    let strand = txs[0].strand;
    let max_exon = txs.iter().map(|t| t.num_exons).max().unwrap();

    let mut collapse_start = Vec::new();
    let mut collapse_end = Vec::new();
    let mut start_wobble = Vec::new();
    let mut end_wobble = Vec::new();
    // support keyed by coord across all priorities
    let mut e_start_support: BTreeMap<i64, IndexSet<String>> = BTreeMap::new();
    let mut e_end_support: BTreeMap<i64, IndexSet<String>> = BTreeMap::new();

    for i in 0..max_exon {
        // priority -> coord -> count
        let mut e_start_dict: BTreeMap<i64, IndexMap<i64, i64>> = BTreeMap::new();
        let mut e_end_dict: BTreeMap<i64, IndexMap<i64, i64>> = BTreeMap::new();

        for t in txs {
            if i >= t.num_exons {
                continue;
            }
            let j = if strand == '+' {
                t.num_exons - 1 - i
            } else {
                i
            };
            let (esp, eep) = exon_priorities(t, i, max_exon, strand);
            let es = t.exon_start_list[j];
            let ee = t.exon_end_list[j];

            // Reset a coordinate's support set when its (priority, coord) pair is
            // first seen at this exon position, mirroring the original's
            // `e_start_trans_dict[e_start] = {}` — a later-priority member (e.g. a
            // no_cap supporter) replaces an earlier-priority one in the support cell.
            let start_bucket = e_start_dict.entry(esp).or_default();
            let start_new = !start_bucket.contains_key(&es);
            *start_bucket.entry(es).or_insert(0) += 1;
            if start_new {
                e_start_support.insert(es, IndexSet::new());
            }
            e_start_support
                .entry(es)
                .or_default()
                .insert(t.uniq_trans_id.clone());

            let end_bucket = e_end_dict.entry(eep).or_default();
            let end_new = !end_bucket.contains_key(&ee);
            *end_bucket.entry(ee).or_insert(0) += 1;
            if end_new {
                e_end_support.insert(ee, IndexSet::new());
            }
            e_end_support
                .entry(ee)
                .or_default()
                .insert(t.uniq_trans_id.clone());
        }

        let bsp = *e_start_dict.keys().next().unwrap();
        let bep = *e_end_dict.keys().next().unwrap();
        let (best_start, long_start, short_start) = best_coord(&e_start_dict[&bsp], true);
        let (best_end, long_end, short_end) = best_coord(&e_end_dict[&bep], false);

        start_wobble.push(short_start - long_start);
        end_wobble.push(long_end - short_end);

        let mut bs = best_start;
        let mut be = best_end;
        if th.longest_ends {
            if i + 1 == max_exon {
                if strand == '+' {
                    bs = long_start;
                } else {
                    be = long_end;
                }
            }
            if i == 0 {
                if strand == '+' {
                    be = long_end;
                } else {
                    bs = long_start;
                }
            }
        }
        collapse_start.push(bs);
        collapse_end.push(be);
    }

    // sort ascending, keep wobble paired
    let mut s: Vec<(i64, i64)> = collapse_start.into_iter().zip(start_wobble).collect();
    s.sort();
    let mut e: Vec<(i64, i64)> = collapse_end.into_iter().zip(end_wobble).collect();
    e.sort();

    let start: Vec<i64> = s.iter().map(|x| x.0).collect();
    let end: Vec<i64> = e.iter().map(|x| x.0).collect();
    let e_start_line: Vec<String> = start
        .iter()
        .map(|c| support_line(&e_start_support, *c))
        .collect();
    let e_end_line: Vec<String> = end
        .iter()
        .map(|c| support_line(&e_end_support, *c))
        .collect();

    Collapsed {
        start,
        end,
        start_wobble: s.iter().map(|x| x.1).collect(),
        end_wobble: e.iter().map(|x| x.1).collect(),
        e_start_support: e_start_line,
        e_end_support: e_end_line,
    }
}

fn support_line(map: &BTreeMap<i64, IndexSet<String>>, coord: i64) -> String {
    let mut v: Vec<String> = map
        .get(&coord)
        .map(|s| s.iter().cloned().collect())
        .unwrap_or_default();
    v.sort();
    v.join(",")
}

/// (best by count with longest tie-break, longest, shortest) for a coord->count map.
/// `is_start`: longest = smallest coord and tie-break to smallest; else longest = largest.
fn best_coord(map: &IndexMap<i64, i64>, is_start: bool) -> (i64, i64, i64) {
    let mut best = if is_start { -1 } else { 0 };
    let mut num = 0i64;
    let mut long = -1i64;
    let mut short = -1i64;
    for (&c, &cnt) in map {
        if cnt > num {
            best = c;
            num = cnt;
        }
        if is_start {
            if long == -1 || c < long {
                long = c;
            }
            if short == -1 || c > short {
                short = c;
            }
        } else {
            if long == -1 || c > long {
                long = c;
            }
            if short == -1 || c < short {
                short = c;
            }
        }
    }
    // tie-break: multiple coords at max count -> pick longest
    let mut num_most = 0;
    let mut most_long = -1i64;
    for (&c, &cnt) in map {
        if cnt == num {
            num_most += 1;
            if is_start {
                if most_long == -1 || most_long > c {
                    most_long = c;
                }
            } else if most_long == -1 || most_long < c {
                most_long = c;
            }
        }
    }
    if num_most > 1 {
        best = most_long;
    }
    (best, long, short)
}

/// Exon start/end priority for exon `i` (3'-indexed). Mirrors the merge
/// `collapse_transcripts` priority assignment for capped reads.
fn exon_priorities(t: &MergeTx, i: usize, max_exon: usize, strand: char) -> (i64, i64) {
    let (sp, jp, ep) = (t.start_priority, t.junction_priority, t.end_priority);
    if i == max_exon - 1 {
        if strand == '+' {
            (sp, jp)
        } else {
            (jp, sp)
        }
    } else if i == t.num_exons - 1 {
        // 5' exon of a shorter read (nocap); forced-last on the 5' side
        if strand == '+' {
            (999, jp)
        } else {
            (jp, 999)
        }
    } else if i > 0 {
        (jp, jp)
    } else if strand == '+' {
        (jp, ep)
    } else {
        (ep, jp)
    }
}

fn rgb_for(members: &[MergeTx], source_colour: &IndexMap<String, usize>) -> String {
    const RGB: [&str; 10] = [
        "255,0,0",
        "255,100,0",
        "255,200,0",
        "200,255,0",
        "0,255,200",
        "0,200,255",
        "0,100,255",
        "0,0,255",
        "100,0,255",
        "200,0,255",
    ];
    let srcs: IndexSet<&str> = members.iter().map(|m| m.source_id.as_str()).collect();
    if srcs.len() > 1 {
        RGB[9].to_string()
    } else {
        let c = source_colour[members[0].source_id.as_str()];
        RGB[c - 1].to_string()
    }
}

fn member_bed_line(t: &MergeTx, final_trans_id: &str) -> String {
    let id_line = format!("{};{}", final_trans_id, t.uniq_trans_id);
    let (mut sizes, mut starts) = (String::new(), String::new());
    for k in 0..t.num_exons {
        if k > 0 {
            sizes.push(',');
            starts.push(',');
        }
        sizes.push_str(&(t.exon_end_list[k] - t.exon_start_list[k]).to_string());
        starts.push_str(&(t.exon_start_list[k] - t.trans_start).to_string());
    }
    [
        t.scaffold.clone(),
        t.trans_start.to_string(),
        t.trans_end.to_string(),
        id_line,
        "40".to_string(),
        t.strand.to_string(),
        t.trans_start.to_string(),
        t.trans_end.to_string(),
        "255,0,0".to_string(),
        t.num_exons.to_string(),
        sizes,
        starts,
    ]
    .join("\t")
}

fn format_trans_report(m: &MergedTrans) -> String {
    let mut sources: Vec<String> = Vec::new();
    for uid in &m.members {
        let s = uid.split('_').next().unwrap_or("").to_string();
        if !sources.contains(&s) {
            sources.push(s);
        }
    }
    sources.sort();
    let join_i = |v: &[i64]| {
        v.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    };
    [
        m.trans_id.clone(),
        m.members.len().to_string(),
        sources.join(","),
        join_i(&m.collapsed.start_wobble),
        join_i(&m.collapsed.end_wobble),
        m.collapsed.e_start_support.join(";"),
        m.collapsed.e_end_support.join(";"),
        m.members.join(","),
    ]
    .join("\t")
}

fn format_gene_report(
    gene_txs: &[MergeTx],
    total_gene_count: i64,
    trans_count: i64,
    src_genes: &[String],
) -> String {
    let gene_id = format!("G{total_gene_count}");
    let mut gene_start = i64::MAX;
    let mut gene_end = 0i64;
    let chrom = gene_txs[0].scaffold.clone();
    let mut src: IndexSet<String> = IndexSet::new();
    for t in gene_txs {
        gene_start = gene_start.min(t.trans_start);
        gene_end = gene_end.max(t.trans_end);
        src.insert(t.source_id.clone());
    }
    let mut src_list: Vec<String> = src.into_iter().collect();
    src_list.sort();

    // source summary: source -> count of source genes
    let mut summary: IndexMap<String, i64> = IndexMap::new();
    for sg in src_genes {
        let s = sg.split('_').next().unwrap_or("").to_string();
        *summary.entry(s).or_insert(0) += 1;
    }
    let summary_line = summary
        .iter()
        .map(|(s, c)| format!("{s}:{c}"))
        .collect::<Vec<_>>()
        .join(",");

    [
        gene_id,
        gene_txs.len().to_string(),
        trans_count.to_string(),
        src_list.join(","),
        chrom,
        gene_start.to_string(),
        gene_end.to_string(),
        src_genes.join(","),
        summary_line,
    ]
    .join("\t")
}

fn sort_merged(merged: &mut [MergedTrans]) {
    merged.sort_by(|a, b| {
        let (ka, kb) = (a.pos_key(), b.pos_key());
        let n = ka.len().max(kb.len());
        for k in 0..n {
            let av = ka.get(k).copied().unwrap_or(0);
            let bv = kb.get(k).copied().unwrap_or(0);
            match av.cmp(&bv) {
                std::cmp::Ordering::Equal => continue,
                o => return o,
            }
        }
        std::cmp::Ordering::Equal
    });
}
