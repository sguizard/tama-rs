//! Per-read error/mismatch computation, ported from `tama_collapse.py`
//! (`mismatch_seq`, `update_variation_dict`, `calc_error_rate`).
//!
//! Walks a CIGAR against the genome to count hard/soft clips, insertions,
//! deletions and mismatches; records variation; and builds the pre/post
//! splice-junction error strings used by the trans report and LDE output.

use std::collections::BTreeMap;
use std::collections::HashMap;

use indexmap::{IndexMap, IndexSet};

use crate::cigar::cigar_list;
use crate::error::Error;

/// Marker appended when no mismatch is found near a splice junction
/// (`no_mismatch_flag` in the original).
pub const NO_MISMATCH_FLAG: &str = "0";

/// IUPAC expansion of a nucleotide code to the concrete bases it can represent.
fn nuc_set(code: u8) -> &'static [u8] {
    match code.to_ascii_uppercase() {
        b'A' => b"A",
        b'T' => b"T",
        b'C' => b"C",
        b'G' => b"G",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"AC",
        b'M' => b"GT",
        b'Y' => b"AG",
        b'R' => b"CT",
        b'V' => b"CGT",
        b'H' => b"AGT",
        b'D' => b"ACT",
        b'B' => b"ACG",
        b'N' => b"ATCG",
        _ => b"",
    }
}

/// Do a genome base and a read base share any concrete nucleotide? Mirrors the
/// set-overlap test in `mismatch_seq`.
fn nuc_match(genome: u8, read: u8) -> bool {
    let g = nuc_set(genome);
    let r = nuc_set(read.to_ascii_uppercase());
    g.iter().any(|c| r.contains(c))
}

/// Result of comparing an aligned genome/query slice.
struct MismatchSlice {
    genome_mismatch: Vec<i64>,
    seq_mismatch: Vec<i64>,
    nuc_mismatch: Vec<String>,
}

/// Compare equal-length genome and query slices, returning per-mismatch genome
/// coordinates, query coordinates, and `query.genome` nucleotide strings.
fn mismatch_seq(
    genome_seq: &[u8],
    query_seq: &[u8],
    genome_pos: i64,
    seq_pos: i64,
) -> Result<MismatchSlice, Error> {
    if genome_seq.len() != query_seq.len() {
        return Err(Error::Invalid(format!(
            "genome slice ({}) != query slice ({})",
            genome_seq.len(),
            query_seq.len()
        )));
    }
    let mut out = MismatchSlice {
        genome_mismatch: Vec::new(),
        seq_mismatch: Vec::new(),
        nuc_mismatch: Vec::new(),
    };
    for i in 0..genome_seq.len() {
        let g = genome_seq[i];
        let r = query_seq[i];
        if !nuc_match(g, r) {
            out.genome_mismatch.push(genome_pos + i as i64);
            out.seq_mismatch.push(seq_pos + i as i64);
            out.nuc_mismatch.push(format!(
                "{}.{}",
                (r as char).to_ascii_uppercase(),
                g as char
            ));
        }
    }
    Ok(out)
}

/// scaffold -> pos -> var_type -> var_seq -> set(read_id). The alt-seq level uses
/// insertion-ordered maps so variant lines emit deterministically.
pub type VariationDict =
    HashMap<String, BTreeMap<i64, IndexMap<char, IndexMap<String, IndexSet<String>>>>>;
/// scaffold -> pos -> set(read_id).
pub type CoverageDict = HashMap<String, BTreeMap<i64, IndexSet<String>>>;

/// Variation and per-position coverage accumulators.
///
/// Insertion order of read IDs is preserved to match the Python dict iteration
/// used when emitting variants.
#[derive(Default)]
pub struct Variation {
    pub dict: VariationDict,
    pub coverage: CoverageDict,
}

impl Variation {
    fn update(&mut self, scaffold: &str, pos: i64, vtype: char, seq: &str, read_id: &str) {
        self.dict
            .entry(scaffold.to_string())
            .or_default()
            .entry(pos)
            .or_default()
            .entry(vtype)
            .or_default()
            .entry(seq.to_string())
            .or_default()
            .insert(read_id.to_string());
        self.coverage
            .entry(scaffold.to_string())
            .or_default()
            .entry(pos)
            .or_default()
            .insert(read_id.to_string());
    }
}

/// Error counts and splice-junction error strings for one read.
#[derive(Debug, Clone, Default)]
pub struct ErrorRate {
    pub h_count: i64,
    pub s_count: i64,
    pub i_count: i64,
    pub d_count: i64,
    pub mis_count: i64,
    pub insertion_list: Vec<i64>,
    pub insertion_length_list: Vec<i64>,
    pub deletion_list: Vec<i64>,
    pub deletion_length_list: Vec<i64>,
    pub mismatch_list: Vec<i64>,
    /// One entry per splice junction (5'->3'), errors within `sj_err_threshold`.
    pub sj_pre_error_list: Vec<String>,
    pub sj_post_error_list: Vec<String>,
}

/// Compute error counts, variation, and splice-junction error strings for a read.
///
/// `genome` is the full (uppercased) scaffold sequence; `seq_list` is the read
/// query sequence. `start_pos` is 1-based (SAM POS). Faithful port of
/// `calc_error_rate`.
#[allow(clippy::too_many_arguments)]
pub fn calc_error_rate(
    start_pos: i64,
    cigar: &str,
    seq_list: &[u8],
    genome: &[u8],
    scaffold: &str,
    read_id: &str,
    sj_err_threshold: i64,
    variation: &mut Variation,
) -> Result<ErrorRate, Error> {
    let ops = cigar_list(cigar)?;
    let lens: Vec<i64> = ops.iter().map(|o| o.len).collect();
    let chars: Vec<u8> = ops.iter().map(|o| o.op).collect();

    let mut r = ErrorRate::default();
    let mut all_nuc_mismatch: Vec<String> = Vec::new();

    let mut genome_pos = start_pos - 1; // 1-based -> 0-based
    let mut seq_pos: i64 = 0;

    let gslice = |a: i64, b: i64| -> &[u8] { &genome[a as usize..b as usize] };
    let qslice = |a: i64, b: i64| -> &[u8] { &seq_list[a as usize..b as usize] };

    for i in 0..ops.len() {
        let flag = chars[i];
        let len = lens[i];
        match flag {
            b'H' => {
                r.h_count += len;
                let var_pos = genome_pos - r.h_count;
                let var_seq = "N".repeat(len as usize);
                variation.update(scaffold, var_pos, 'H', &var_seq, read_id);
            }
            b'S' => {
                r.s_count += len;
                let seq_start = seq_pos;
                let seq_end = seq_pos + len;
                let mut var_pos = genome_pos - r.s_count;
                let var_seq = String::from_utf8_lossy(qslice(seq_start, seq_end)).into_owned();
                if i == ops.len() - 1 {
                    var_pos = genome_pos + 1 - r.s_count;
                }
                variation.update(scaffold, var_pos, 'S', &var_seq, read_id);
                seq_pos += len;
            }
            b'M' => {
                let seq_start = seq_pos;
                let seq_end = seq_pos + len;
                let genome_start = genome_pos;
                let genome_end = genome_pos + len;
                let ms = mismatch_seq(
                    gslice(genome_start, genome_end),
                    qslice(seq_start, seq_end),
                    genome_start,
                    seq_start,
                )?;
                r.mis_count += ms.genome_mismatch.len() as i64;
                for k in 0..ms.seq_mismatch.len() {
                    let var_pos = ms.genome_mismatch[k];
                    let seq_var_pos = ms.seq_mismatch[k];
                    let var_seq = (seq_list[seq_var_pos as usize] as char).to_string();
                    variation.update(scaffold, var_pos, 'M', &var_seq, read_id);
                }
                r.mismatch_list.extend(&ms.genome_mismatch);
                all_nuc_mismatch.extend(ms.nuc_mismatch);
                seq_pos += len;
                genome_pos += len;
            }
            b'I' => {
                let seq_start = seq_pos;
                let seq_end = seq_pos + len;
                let var_seq = String::from_utf8_lossy(qslice(seq_start, seq_end)).into_owned();
                variation.update(scaffold, genome_pos, 'I', &var_seq, read_id);
                seq_pos += len;
                r.i_count += len;
                r.insertion_list.push(genome_pos);
                r.insertion_length_list.push(len);
            }
            b'D' => {
                r.deletion_list.push(genome_pos);
                r.deletion_length_list.push(len);
                variation.update(scaffold, genome_pos, 'D', &len.to_string(), read_id);
                genome_pos += len;
                r.d_count += len;
            }
            b'N' => {
                let (pre, post) = sj_errors(
                    &lens,
                    &chars,
                    i,
                    genome_pos,
                    seq_pos,
                    &r.mismatch_list,
                    &all_nuc_mismatch,
                    genome,
                    seq_list,
                    scaffold,
                    sj_err_threshold,
                );
                r.sj_pre_error_list.push(pre);
                r.sj_post_error_list.push(post);
                genome_pos += len; // intron
            }
            other => {
                return Err(Error::Cigar(format!(
                    "unexpected CIGAR op {:?} in {cigar:?}",
                    other as char
                )));
            }
        }
    }

    r.h_count = r.h_count.max(0);
    Ok(r)
}

/// Build the pre/post splice-junction error strings for the `N` op at index `i`.
/// Faithful port of the two `while` loops in `calc_error_rate`.
#[allow(clippy::too_many_arguments)]
fn sj_errors(
    lens: &[i64],
    chars: &[u8],
    i: usize,
    genome_pos: i64,
    seq_pos: i64,
    all_genome_mismatch: &[i64],
    all_nuc_mismatch: &[String],
    genome: &[u8],
    seq_list: &[u8],
    _scaffold: &str,
    sj_err_threshold: i64,
) -> (String, String) {
    // ---- pre (before the junction), walked backwards ----
    let mut pre: Vec<String> = Vec::new();
    {
        let mut prev_cig_flag = chars[i - 1];
        let mut prev_cig_length = lens[i - 1];
        let mut this_mismatch_length = 0i64;
        let mut prev_total = 0i64;
        let mut mm_index: i64 = all_genome_mismatch.len() as i64 - 1;
        let mut prev_cig_index: i64 = i as i64 - 1;
        let mut this_cig_genome_pos = genome_pos;
        let mut prev_sj_flag = false;

        while prev_total <= sj_err_threshold && !prev_sj_flag {
            this_mismatch_length += prev_cig_length;
            prev_total = this_mismatch_length;

            if prev_cig_flag == b'M' {
                if all_genome_mismatch.is_empty() {
                    if prev_total <= sj_err_threshold {
                        pre.push(format!("{}{}", prev_cig_length, prev_cig_flag as char));
                    }
                } else if mm_index >= 0 {
                    let mut last_gm = all_genome_mismatch[mm_index as usize];
                    let mut last_nm = all_nuc_mismatch[mm_index as usize].clone();
                    let mut dist = genome_pos - last_gm;
                    let mut add_count = 0;
                    while last_gm >= this_cig_genome_pos - prev_cig_length
                        && last_gm <= this_cig_genome_pos
                        && dist <= sj_err_threshold
                    {
                        pre.push(format!("{}.{}", dist, last_nm));
                        add_count += 1;
                        mm_index -= 1;
                        if mm_index < 0 {
                            break;
                        }
                        last_gm = all_genome_mismatch[mm_index as usize];
                        last_nm = all_nuc_mismatch[mm_index as usize].clone();
                        dist = genome_pos - last_gm;
                    }
                    if add_count == 0 && prev_total <= sj_err_threshold {
                        pre.push(format!("{}{}", prev_cig_length, prev_cig_flag as char));
                    }
                } else if prev_total <= sj_err_threshold {
                    pre.push(format!("{}{}", prev_cig_length, prev_cig_flag as char));
                }
            } else if prev_cig_flag != b'N' {
                prev_cig_flag = chars[prev_cig_index as usize];
                prev_cig_length = lens[prev_cig_index as usize];
                pre.push(format!("{}{}", prev_cig_length, prev_cig_flag as char));
            } else {
                prev_sj_flag = true;
            }

            prev_cig_index -= 1;
            if prev_cig_index < 0 {
                break;
            }
            this_cig_genome_pos -= prev_cig_length;
            prev_cig_flag = chars[prev_cig_index as usize];
            prev_cig_length = lens[prev_cig_index as usize];
        }
        if pre.is_empty() {
            pre.push(NO_MISMATCH_FLAG.to_string());
        }
    }

    // ---- post (after the junction), walked forwards ----
    let mut post: Vec<String> = Vec::new();
    {
        let post_genome_pos = genome_pos + lens[i]; // after intron
        let mut next_cig_flag = chars[i + 1];
        let mut next_cig_length = lens[i + 1];
        let mut this_mismatch_length = 0i64;
        let mut next_total = 0i64;
        let mut next_cig_index = i + 1;
        let mut this_next_seq_pos = seq_pos;
        let mut this_genome_pos = post_genome_pos;
        let mut next_sj_flag = false;

        while next_total <= sj_err_threshold && !next_sj_flag {
            this_mismatch_length += next_cig_length;
            next_total = this_mismatch_length;

            if next_cig_flag == b'M' {
                let seq_start = this_next_seq_pos;
                let seq_end = this_next_seq_pos + next_cig_length;
                let genome_start = this_genome_pos;
                let genome_end = this_genome_pos + next_cig_length;
                let ms = mismatch_seq(
                    &genome[genome_start as usize..genome_end as usize],
                    &seq_list[seq_start as usize..seq_end as usize],
                    genome_start,
                    seq_start,
                )
                .expect("equal-length slices");

                if ms.genome_mismatch.is_empty() {
                    if next_total <= sj_err_threshold {
                        post.push(format!("{}{}", next_cig_length, next_cig_flag as char));
                    }
                } else {
                    let mut gm_index = 0usize;
                    let mut dist = ms.genome_mismatch[0] - post_genome_pos;
                    let mut add_count = 0;
                    while dist <= this_mismatch_length && dist < sj_err_threshold {
                        let next_gm = ms.genome_mismatch[gm_index];
                        let next_nm = &ms.nuc_mismatch[gm_index];
                        dist = next_gm - post_genome_pos;
                        if dist <= sj_err_threshold {
                            post.push(format!("{}.{}", dist, next_nm));
                            add_count += 1;
                        }
                        gm_index += 1;
                        if gm_index >= ms.genome_mismatch.len() {
                            break;
                        }
                    }
                    if add_count == 0 && next_total <= sj_err_threshold {
                        post.push(format!("{}{}", next_cig_length, next_cig_flag as char));
                    }
                }
            } else if next_cig_flag != b'N' {
                next_cig_flag = chars[next_cig_index];
                next_cig_length = lens[next_cig_index];
                post.push(format!("{}{}", next_cig_length, next_cig_flag as char));
            } else {
                next_sj_flag = true;
            }

            if next_cig_flag != b'D' {
                this_next_seq_pos += next_cig_length;
            }
            if next_cig_flag != b'I' {
                this_genome_pos += next_cig_length;
            }
            next_cig_index += 1;
            if next_cig_index >= chars.len() {
                break;
            }
            next_cig_flag = chars[next_cig_index];
            next_cig_length = lens[next_cig_index];
        }
        if post.is_empty() {
            post.push(NO_MISMATCH_FLAG.to_string());
        }
    }

    pre.reverse();
    (pre.join("_"), post.join("_"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_simple_perfect_match() {
        // genome = query, single 5M exon, no errors.
        let genome = b"ACGTACGTACGT";
        let query = b"ACGT";
        let mut var = Variation::default();
        let r =
            calc_error_rate(1, "4M", query, genome, "chr", "r1", 10, &mut var).unwrap();
        assert_eq!((r.h_count, r.s_count, r.i_count, r.d_count, r.mis_count), (0, 0, 0, 0, 0));
    }

    #[test]
    fn counts_mismatch_and_indels() {
        //          pos1..: A C G T
        let genome = b"ACGTACGT";
        // 2M (AC), 1I (T), 2M vs genome GT but read has GA -> 1 mismatch at last base
        let query = b"ACTGA";
        let mut var = Variation::default();
        // CIGAR: 2M1I2M  query: AC | T | GA ; genome: AC .. GT
        let r =
            calc_error_rate(1, "2M1I2M", query, genome, "chr", "r1", 10, &mut var).unwrap();
        assert_eq!(r.i_count, 1);
        assert_eq!(r.mis_count, 1); // A vs T at genome pos 4 (0-based 3)
    }

    #[test]
    fn soft_and_hard_clips_counted() {
        let genome = b"ACGTACGTACGT";
        let query = b"NNACGT"; // 2S then 4M
        let mut var = Variation::default();
        let r =
            calc_error_rate(1, "2S4M", query, genome, "chr", "r1", 10, &mut var).unwrap();
        assert_eq!(r.s_count, 2);
        assert_eq!(r.mis_count, 0);
    }
}
