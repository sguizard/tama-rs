//! Transcript collapsing, ported from `tama_collapse.py`.
//!
//! Covers the post-SAM processing: comparing transcripts (`compare_transcripts`),
//! grouping collapsible transcripts within a gene (`simplify_gene_capped`, as
//! connected components over the `same_transcript` relation), merging a group's
//! exon coordinates (`collapse_transcripts`), and sorting/deduplicating the final
//! models (`sort_transcripts`). Splice-junction error priority
//! (`sj_error_priority_finder`) is included so the trans report matches.
//!
//! The no-cap grouping path (`simplify_gene_nocap`) is not yet ported; callers
//! should reject `no_cap` until it lands.

use std::collections::BTreeMap;

use indexmap::{IndexMap, IndexSet};

pub const NO_MISMATCH_FLAG: &str = "0";

/// End-collapse behaviour (`-e`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ends {
    CommonEnds,
    LongestEnds,
}

/// 5' cap handling (`-x`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cap {
    Capped,
    NoCap,
}

/// Thresholds and flags controlling collapsing.
#[derive(Debug, Clone, Copy)]
pub struct CollapseParams {
    pub fiveprime_threshold: i64,
    pub threeprime_threshold: i64,
    pub exon_diff_threshold: i64,
    pub cap: Cap,
    pub ends: Ends,
    pub sj_priority: bool,
    pub merge_dup: bool,
}

/// A read-level transcript model used during collapsing.
#[derive(Debug, Clone)]
pub struct Transcript {
    pub cluster_id: String,
    pub scaff_name: String,
    pub strand: char,
    pub start_pos: i64,
    pub end_pos: i64,
    pub exon_start_list: Vec<i64>,
    pub exon_end_list: Vec<i64>,
    pub sj_pre_error_list: Vec<String>,
    pub sj_post_error_list: Vec<String>,
    // carried for reporting only
    pub percent_coverage: f64,
    pub percent_identity: f64,
}

impl Transcript {
    pub fn num_exons(&self) -> usize {
        self.exon_start_list.len()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MatchFlag {
    Perfect,
    Wobbly,
    NoMatch,
}

fn fuzzy_match(c1: i64, c2: i64, threshold: i64) -> (MatchFlag, i64) {
    if c1 == c2 {
        (MatchFlag::Perfect, 0)
    } else {
        let d = c1 - c2;
        if d.abs() <= threshold {
            (MatchFlag::Wobbly, d)
        } else {
            (MatchFlag::NoMatch, d)
        }
    }
}

/// Result of comparing two transcripts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompFlag {
    None,
    DiffTranscripts,
    SameTranscript,
    SameThreePrimeSameExons,
    SameThreePrimeDiffExons,
}

/// Compare two transcripts for collapsibility. Faithful port of
/// `compare_transcripts`; returns the comparison flag (the only field the
/// grouping needs).
// The no-cap branches keep the original's separate a/b conditions even though
// some set the same flag, to stay aligned with the Python source.
#[allow(clippy::if_same_then_else)]
pub fn compare_transcripts(
    a: &Transcript,
    b: &Transcript,
    p: &CollapseParams,
    strand: char,
) -> CompFlag {
    let a_id = &a.cluster_id;
    let b_id = &b.cluster_id;
    let e_start = &a.exon_start_list;
    let o_e_start = &b.exon_start_list;
    let e_end = &a.exon_end_list;
    let o_e_end = &b.exon_end_list;

    let diff_num_exon_flag: i64 = if e_start.len() != o_e_start.len() {
        e_start.len() as i64 - o_e_start.len() as i64
    } else {
        0
    };

    let capped = p.cap == Cap::Capped;
    if capped && diff_num_exon_flag != 0 {
        return CompFlag::DiffTranscripts;
    }

    let min_exon_num = e_start.len().min(o_e_start.len());
    // determine short/long for equal-length case (nocap semantics)
    let mut short_trans: String;
    if e_start.len() == o_e_start.len() {
        if strand == '+' {
            if e_start[0] < o_e_start[0] {
                short_trans = b_id.clone();
            } else if e_start[0] > o_e_start[0] {
                short_trans = a_id.clone();
            } else {
                short_trans = "same".to_string();
            }
        } else if e_end[0] > o_e_end[0] {
            short_trans = b_id.clone();
        } else if e_end[0] < o_e_end[0] {
            short_trans = a_id.clone();
        } else {
            short_trans = "same".to_string();
        }
    } else if e_start.len() > o_e_start.len() {
        short_trans = b_id.clone();
    } else {
        short_trans = a_id.clone();
    }

    let mut flag = CompFlag::None;
    let mut all_match = true;

    for i in 0..min_exon_num {
        // Python indexes each list from its own 3' end with `-(i+1)` on the plus
        // strand, or from the 5' end with `i` on the minus strand.
        let (ja, jb) = if strand == '+' {
            (e_start.len() - 1 - i, o_e_start.len() - 1 - i)
        } else {
            (i, i)
        };
        let es = e_start[ja];
        let oes = o_e_start[jb];
        let ee = e_end[ja];
        let oee = o_e_end[jb];

        // micro-exon: non-overlap -> different transcripts
        if es >= oee || oes >= ee {
            flag = CompFlag::DiffTranscripts;
            continue;
        }

        let mut start_threshold = p.exon_diff_threshold;
        let mut end_threshold = p.exon_diff_threshold;
        if strand == '+' {
            if i == 0 {
                end_threshold = p.threeprime_threshold;
            }
            if diff_num_exon_flag == 0 && i == min_exon_num - 1 {
                start_threshold = p.fiveprime_threshold;
            }
        } else {
            if i == 0 {
                start_threshold = p.threeprime_threshold;
            }
            if diff_num_exon_flag == 0 && i == min_exon_num - 1 {
                end_threshold = p.fiveprime_threshold;
            }
        }

        let (start_match, start_diff) = fuzzy_match(es, oes, start_threshold);
        let (end_match, end_diff) = fuzzy_match(ee, oee, end_threshold);

        let nocap_last = p.cap == Cap::NoCap && i == min_exon_num - 1;
        if nocap_last {
            if strand == '+' {
                if end_match == MatchFlag::NoMatch {
                    all_match = false;
                } else if start_match == MatchFlag::NoMatch && all_match {
                    if *b_id == short_trans && start_diff < 0 && diff_num_exon_flag != 0 {
                        flag = CompFlag::SameThreePrimeDiffExons;
                    } else if *a_id == short_trans && start_diff > 0 && diff_num_exon_flag != 0 {
                        flag = CompFlag::SameThreePrimeDiffExons;
                    } else if start_diff < 0 && diff_num_exon_flag == 0 {
                        short_trans = b_id.clone();
                        flag = CompFlag::SameThreePrimeSameExons;
                    } else if start_diff > 0 && diff_num_exon_flag == 0 {
                        short_trans = a_id.clone();
                        flag = CompFlag::SameThreePrimeSameExons;
                    } else {
                        all_match = false;
                    }
                }
            } else if start_match == MatchFlag::NoMatch {
                all_match = false;
            } else if end_match == MatchFlag::NoMatch && all_match {
                if *b_id == short_trans && end_diff > 0 && diff_num_exon_flag != 0 {
                    flag = CompFlag::SameThreePrimeDiffExons;
                } else if *a_id == short_trans && end_diff < 0 && diff_num_exon_flag != 0 {
                    flag = CompFlag::SameThreePrimeDiffExons;
                } else if end_diff > 0 && diff_num_exon_flag == 0 {
                    short_trans = b_id.clone();
                    flag = CompFlag::SameThreePrimeSameExons;
                } else if end_diff < 0 && diff_num_exon_flag == 0 {
                    short_trans = a_id.clone();
                    flag = CompFlag::SameThreePrimeSameExons;
                } else {
                    all_match = false;
                }
            }
        } else {
            if start_match == MatchFlag::NoMatch {
                all_match = false;
            }
            if end_match == MatchFlag::NoMatch {
                all_match = false;
            }
        }
    }

    if flag == CompFlag::None {
        if all_match {
            if diff_num_exon_flag == 0 && capped {
                flag = CompFlag::SameTranscript;
            } else if diff_num_exon_flag == 0 && !capped {
                flag = CompFlag::SameThreePrimeSameExons;
            } else {
                flag = CompFlag::SameThreePrimeDiffExons;
            }
        } else {
            flag = CompFlag::DiffTranscripts;
        }
    }

    flag
}

/// Splice-junction error priority + error strings for exon `i` (0 = 3'-most).
/// Faithful port of `sj_error_priority_finder`.
pub fn sj_error_priority_finder(
    t: &Transcript,
    i: usize,
    max_exon_num: usize,
    sj_priority: bool,
) -> (i64, i64, String, String) {
    let exon_num = t.exon_start_list.len();
    let pre = &t.sj_pre_error_list;
    let post = &t.sj_post_error_list;
    let mut e_start_priority;
    let mut e_end_priority;
    let mut e_start_err = "na".to_string();
    let mut e_end_err = "na".to_string();
    let delim = ">";

    let start_pri = |a: &str, b: &str| -> i64 {
        match (b == NO_MISMATCH_FLAG, a == NO_MISMATCH_FLAG) {
            (true, true) => 0,
            (true, false) => 1,
            (false, true) => 2,
            (false, false) => 3,
        }
    };
    let end_pri = |pre_e: &str, post_e: &str| -> i64 {
        let pn = pre_e == NO_MISMATCH_FLAG;
        let qn = post_e == NO_MISMATCH_FLAG;
        match (qn, pn) {
            (true, true) => 0,
            (false, true) => 1,
            (true, false) => 2,
            (false, false) => 3,
        }
    };

    if exon_num > 1 {
        if t.strand == '+' {
            if i == 0 {
                e_end_priority = 0;
                let sj_pre = &pre[pre.len() - 1];
                let sj_post = &post[post.len() - 1];
                e_start_priority = start_pri(sj_pre, sj_post);
                e_start_err = format!("{}{}{}", sj_pre, delim, sj_post);
            } else if i < max_exon_num - 1 {
                let sj_pre_s = &pre[pre.len() - 1 - i];
                let sj_post_s = &post[post.len() - 1 - i];
                let sj_pre_e = &pre[pre.len() - i];
                let sj_post_e = &post[post.len() - i];
                e_start_priority = start_pri(sj_pre_s, sj_post_s);
                e_end_priority = end_pri(sj_pre_e, sj_post_e);
                e_start_err = format!("{}{}{}", sj_pre_s, delim, sj_post_s);
                e_end_err = format!("{}{}{}", sj_pre_e, delim, sj_post_e);
            } else {
                e_start_priority = 0;
                let sj_pre_e = &pre[pre.len() - i];
                let sj_post_e = &post[post.len() - i];
                e_end_priority = end_pri(sj_pre_e, sj_post_e);
                e_end_err = format!("{}{}{}", sj_pre_e, delim, sj_post_e);
            }
        } else if i == 0 {
            e_start_priority = 0;
            let sj_post = &post[0];
            let sj_pre = &pre[0];
            e_end_priority = end_pri(sj_pre, sj_post);
            e_end_err = format!("{}{}{}", sj_pre, delim, sj_post);
        } else if i < max_exon_num - 1 {
            let sj_pre_s = &pre[i - 1];
            let sj_post_s = &post[i - 1];
            let sj_pre_e = &pre[i];
            let sj_post_e = &post[i];
            e_start_priority = start_pri(sj_pre_s, sj_post_s);
            e_end_priority = end_pri(sj_pre_e, sj_post_e);
            e_start_err = format!("{}{}{}", sj_pre_s, delim, sj_post_s);
            e_end_err = format!("{}{}{}", sj_pre_e, delim, sj_post_e);
        } else {
            e_end_priority = 0;
            let sj_pre_s = &pre[i - 1];
            let sj_post_s = &post[i - 1];
            e_start_priority = start_pri(sj_pre_s, sj_post_s);
            e_start_err = format!("{}{}{}", sj_pre_s, delim, sj_post_s);
        }
    } else {
        e_start_priority = 0;
        e_end_priority = 0;
    }

    if !sj_priority {
        e_start_priority = 0;
        e_end_priority = 0;
    }

    (e_start_priority, e_end_priority, e_start_err, e_end_err)
}

/// Collapsed exon coordinates and diagnostics for a group of transcripts.
pub struct CollapseResult {
    pub collapse_start_list: Vec<i64>,
    pub collapse_end_list: Vec<i64>,
    pub start_wobble_list: Vec<i64>,
    pub end_wobble_list: Vec<i64>,
    pub collapse_sj_start_err_list: Vec<i64>,
    pub collapse_sj_end_err_list: Vec<i64>,
    pub collapse_start_error_nuc_list: Vec<String>,
    pub collapse_end_error_nuc_list: Vec<String>,
}

/// Merge a group of collapsible transcripts into one model. Faithful port of
/// `collapse_transcripts`.
pub fn collapse_transcripts(trans: &[Transcript], p: &CollapseParams) -> CollapseResult {
    let strand = trans[0].strand;
    let num_trans = trans.len();
    let max_exon_num = trans.iter().map(|t| t.num_exons()).max().unwrap_or(0);

    let mut collapse_start_list = Vec::new();
    let mut collapse_end_list = Vec::new();
    let mut start_wobble_list = Vec::new();
    let mut end_wobble_list = Vec::new();
    let mut collapse_sj_start_err_list = Vec::new();
    let mut collapse_sj_end_err_list = Vec::new();
    let mut collapse_start_error_nuc_list = Vec::new();
    let mut collapse_end_error_nuc_list = Vec::new();

    for i in 0..max_exon_num {
        // priority -> value -> count
        let mut e_start_dict: BTreeMap<i64, IndexMap<i64, i64>> = BTreeMap::new();
        let mut e_end_dict: BTreeMap<i64, IndexMap<i64, i64>> = BTreeMap::new();
        let mut err_start: BTreeMap<i64, IndexMap<i64, IndexSet<String>>> = BTreeMap::new();
        let mut err_end: BTreeMap<i64, IndexMap<i64, IndexSet<String>>> = BTreeMap::new();
        let mut e_start_range: Vec<i64> = Vec::new();
        let mut e_end_range: Vec<i64> = Vec::new();

        for t in trans {
            let this_max = t.num_exons();
            if i >= this_max {
                continue;
            }
            let j = if strand == '+' { this_max - 1 - i } else { i };
            let (sp, ep, se, ee_err) = sj_error_priority_finder(t, i, this_max, p.sj_priority);

            let mut e_start = t.exon_start_list[j];
            let mut e_end = t.exon_end_list[j];

            if p.cap == Cap::NoCap && i == this_max - 1 && i < max_exon_num - 1 {
                if strand == '+' {
                    e_start = -1;
                } else {
                    e_end = -1;
                }
            }

            if e_start != -1 {
                *e_start_dict.entry(sp).or_default().entry(e_start).or_insert(0) += 1;
                err_start
                    .entry(sp)
                    .or_default()
                    .entry(e_start)
                    .or_default()
                    .insert(se.clone());
                e_start_range.push(e_start);
            }
            if e_end != -1 {
                *e_end_dict.entry(ep).or_default().entry(e_end).or_insert(0) += 1;
                err_end
                    .entry(ep)
                    .or_default()
                    .entry(e_end)
                    .or_default()
                    .insert(ee_err.clone());
                e_end_range.push(e_end);
            }
        }

        let best_start_priority = *e_start_dict.keys().next().unwrap();
        let best_end_priority = *e_end_dict.keys().next().unwrap();
        collapse_sj_start_err_list.push(best_start_priority);
        collapse_sj_end_err_list.push(best_end_priority);

        let start_map = &e_start_dict[&best_start_priority];
        let end_map = &e_end_dict[&best_end_priority];

        // best start: highest count, ties -> smallest coord (longest)
        let mut long_e_start = -1i64;
        let mut best_e_start = -1i64;
        let mut num_starts = 0i64;
        for (&s, &c) in start_map {
            if c > num_starts {
                best_e_start = s;
                num_starts = c;
            }
            if long_e_start == -1 || s < long_e_start {
                long_e_start = s;
            }
        }
        let mut most_long_e_start = -1i64;
        let mut num_most_starts = 0;
        for (&s, &c) in start_map {
            if c == num_starts {
                num_most_starts += 1;
                if most_long_e_start == -1 || most_long_e_start > s {
                    most_long_e_start = s;
                }
            }
        }
        if num_most_starts > 1 {
            best_e_start = most_long_e_start;
        }

        e_start_range.sort_unstable();
        start_wobble_list.push(e_start_range[e_start_range.len() - 1] - e_start_range[0]);

        // best end: highest count, ties -> largest coord (longest)
        let mut long_e_end = -1i64;
        let mut best_e_end = 0i64;
        let mut num_ends = 0i64;
        for (&e, &c) in end_map {
            if c > num_ends {
                best_e_end = e;
                num_ends = c;
            }
            if long_e_end == -1 || e > long_e_end {
                long_e_end = e;
            }
        }
        let mut most_long_e_end = -1i64;
        let mut num_most_ends = 0;
        for (&e, &c) in end_map {
            if c == num_ends {
                num_most_ends += 1;
                if most_long_e_end == -1 || most_long_e_end < e {
                    most_long_e_end = e;
                }
            }
        }
        if num_most_ends > 1 {
            best_e_end = most_long_e_end;
        }

        e_end_range.sort_unstable();
        end_wobble_list.push(e_end_range[e_end_range.len() - 1] - e_end_range[0]);

        if p.cap == Cap::NoCap && i + 1 == max_exon_num {
            if strand == '+' {
                best_e_start = long_e_start;
            } else {
                best_e_end = long_e_end;
            }
        }
        if p.ends == Ends::LongestEnds {
            if i + 1 == max_exon_num {
                if strand == '+' {
                    best_e_start = long_e_start;
                } else {
                    best_e_end = long_e_end;
                }
            }
            if i == 0 {
                if strand == '+' {
                    best_e_end = long_e_end;
                } else {
                    best_e_start = long_e_start;
                }
            }
        }

        let _ = num_trans;
        collapse_start_list.push(best_e_start);
        collapse_end_list.push(best_e_end);

        let se_line = err_start[&best_start_priority][&best_e_start]
            .iter()
            .cloned()
            .collect::<Vec<_>>()
            .join("-");
        let ee_line = err_end[&best_end_priority][&best_e_end]
            .iter()
            .cloned()
            .collect::<Vec<_>>()
            .join("-");
        collapse_start_error_nuc_list.push(se_line);
        collapse_end_error_nuc_list.push(ee_line);
    }

    // sort coords ascending, keeping wobble paired
    let mut starts: Vec<(i64, i64)> = collapse_start_list
        .iter()
        .cloned()
        .zip(start_wobble_list.iter().cloned())
        .collect();
    starts.sort();
    let mut ends: Vec<(i64, i64)> = collapse_end_list
        .iter()
        .cloned()
        .zip(end_wobble_list.iter().cloned())
        .collect();
    ends.sort();

    CollapseResult {
        collapse_start_list: starts.iter().map(|x| x.0).collect(),
        collapse_end_list: ends.iter().map(|x| x.0).collect(),
        start_wobble_list: starts.iter().map(|x| x.1).collect(),
        end_wobble_list: ends.iter().map(|x| x.1).collect(),
        collapse_sj_start_err_list,
        collapse_sj_end_err_list,
        collapse_start_error_nuc_list,
        collapse_end_error_nuc_list,
    }
}

/// Group a gene's transcripts into collapsible sets (capped mode). Returns the
/// index groups. Equivalent to `simplify_gene_capped`: connected components over
/// the `same_transcript` relation.
pub fn simplify_gene_capped(trans: &[Transcript], p: &CollapseParams) -> Vec<Vec<usize>> {
    let n = trans.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while parent[r] != r {
            r = parent[r];
        }
        let mut c = x;
        while parent[c] != r {
            let nx = parent[c];
            parent[c] = r;
            c = nx;
        }
        r
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if compare_transcripts(&trans[i], &trans[j], p, trans[i].strand)
                == CompFlag::SameTranscript
            {
                let ri = find(&mut parent, i);
                let rj = find(&mut parent, j);
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }
    let mut groups: IndexMap<usize, Vec<usize>> = IndexMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    groups.into_values().collect()
}

/// A collapsed (merged) transcript ready for output.
#[derive(Debug, Clone)]
pub struct Merged {
    pub trans_id: String,
    pub scaff_name: String,
    pub strand: char,
    pub start_pos: i64,
    pub end_pos: i64,
    pub num_exons: usize,
    pub collapse_start_list: Vec<i64>,
    pub collapse_end_list: Vec<i64>,
    pub start_wobble_list: Vec<i64>,
    pub end_wobble_list: Vec<i64>,
    pub collapse_sj_start_err_list: Vec<i64>,
    pub collapse_sj_end_err_list: Vec<i64>,
    pub collapse_start_error_nuc_list: Vec<String>,
    pub collapse_end_error_nuc_list: Vec<String>,
    /// Member read cluster IDs, in insertion order.
    pub merged_trans: Vec<String>,
}

impl Merged {
    /// TAMA BED12 line for this merged model (mirrors `Merged.format_bed_line`).
    pub fn format_bed_line(&self) -> String {
        let gene_id = self.trans_id.split('.').next().unwrap_or(&self.trans_id);
        let id_line = format!("{};{}", gene_id, self.trans_id);
        let mut sizes = String::new();
        let mut starts = String::new();
        for k in 0..self.num_exons {
            if k > 0 {
                sizes.push(',');
                starts.push(',');
            }
            let es = self.collapse_start_list[k];
            let ee = self.collapse_end_list[k];
            sizes.push_str(&(ee - es).to_string());
            starts.push_str(&(es - self.start_pos).to_string());
        }
        [
            self.scaff_name.clone(),
            (self.start_pos - 1).to_string(),
            (self.end_pos - 1).to_string(),
            id_line,
            "40".to_string(),
            self.strand.to_string(),
            (self.start_pos - 1).to_string(),
            (self.end_pos - 1).to_string(),
            "255,0,0".to_string(),
            self.num_exons.to_string(),
            sizes,
            starts,
        ]
        .join("\t")
    }

    /// Position key used for sorting/dedup: `[start, end, e_start0, e_end0, ...]`.
    fn pos_key(&self) -> Vec<i64> {
        let mut starts = self.collapse_start_list.clone();
        let mut ends = self.collapse_end_list.clone();
        starts.sort_unstable();
        ends.sort_unstable();
        let mut key = vec![starts[0], ends[ends.len() - 1]];
        for k in 0..starts.len() {
            key.push(starts[k]);
            key.push(ends[k]);
        }
        key
    }
}

/// Sort merged transcripts by genomic position (start, end, exon coords),
/// deduplicating identical models. Faithful to `sort_transcripts` +
/// `sort_pos_trans_list` (padded lexicographic order). Duplicate models are
/// merged when `merge_dup` is set.
pub fn sort_transcripts(mut merged: Vec<Merged>, p: &CollapseParams) -> Vec<Merged> {
    if p.merge_dup {
        let mut seen: IndexMap<Vec<i64>, usize> = IndexMap::new();
        let mut out: Vec<Merged> = Vec::new();
        for m in merged.drain(..) {
            let key = m.pos_key();
            if let Some(&idx) = seen.get(&key) {
                let existing = &mut out[idx];
                for c in m.merged_trans {
                    if !existing.merged_trans.contains(&c) {
                        existing.merged_trans.push(c);
                    }
                }
            } else {
                seen.insert(key, out.len());
                out.push(m);
            }
        }
        merged = out;
    }
    merged.sort_by(|a, b| pad_cmp(&a.pos_key(), &b.pos_key()));
    merged
}

/// Compare two position keys with 0-padding to equal length (lexicographic).
fn pad_cmp(a: &[i64], b: &[i64]) -> std::cmp::Ordering {
    let n = a.len().max(b.len());
    for k in 0..n {
        let av = a.get(k).copied().unwrap_or(0);
        let bv = b.get(k).copied().unwrap_or(0);
        match av.cmp(&bv) {
            std::cmp::Ordering::Equal => continue,
            other => return other,
        }
    }
    std::cmp::Ordering::Equal
}

#[cfg(test)]
mod tests {
    use super::*;

    fn t(id: &str, strand: char, starts: &[i64], ends: &[i64]) -> Transcript {
        let n = starts.len();
        Transcript {
            cluster_id: id.to_string(),
            scaff_name: "chr".to_string(),
            strand,
            start_pos: starts[0],
            end_pos: ends[n - 1],
            exon_start_list: starts.to_vec(),
            exon_end_list: ends.to_vec(),
            sj_pre_error_list: vec![NO_MISMATCH_FLAG.to_string(); n.saturating_sub(1)],
            sj_post_error_list: vec![NO_MISMATCH_FLAG.to_string(); n.saturating_sub(1)],
            percent_coverage: 100.0,
            percent_identity: 100.0,
        }
    }

    fn params() -> CollapseParams {
        CollapseParams {
            fiveprime_threshold: 10,
            threeprime_threshold: 10,
            exon_diff_threshold: 10,
            cap: Cap::Capped,
            ends: Ends::CommonEnds,
            sj_priority: false,
            merge_dup: true,
        }
    }

    #[test]
    fn identical_transcripts_are_same() {
        let a = t("a", '+', &[100, 500], &[200, 600]);
        let b = t("b", '+', &[103, 500], &[200, 597]);
        assert_eq!(
            compare_transcripts(&a, &b, &params(), '+'),
            CompFlag::SameTranscript
        );
    }

    #[test]
    fn different_exon_count_capped_differs() {
        let a = t("a", '+', &[100, 500, 900], &[200, 600, 1000]);
        let b = t("b", '+', &[100, 500], &[200, 600]);
        assert_eq!(
            compare_transcripts(&a, &b, &params(), '+'),
            CompFlag::DiffTranscripts
        );
    }

    #[test]
    fn collapse_picks_common_ends() {
        // three transcripts, common internal SJ, wobbling 3' end.
        let a = t("a", '+', &[100, 500], &[205, 600]);
        let b = t("b", '+', &[100, 500], &[200, 600]);
        let c = t("c", '+', &[100, 500], &[200, 600]);
        let group = simplify_gene_capped(&[a.clone(), b.clone(), c.clone()], &params());
        assert_eq!(group.len(), 1);
        let all = [a, b, c];
        let r = collapse_transcripts(&all, &params());
        // exon 0 (3') end: 200 appears twice, 205 once -> 200 chosen
        assert_eq!(r.collapse_start_list, vec![100, 500]);
        assert_eq!(r.collapse_end_list, vec![200, 600]);
    }
}
