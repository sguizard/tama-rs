//! `tama cleanup` — sequence cleanup tools.

use std::collections::BTreeMap;
use std::io::{BufRead, Write};

use clap::{Args as ClapArgs, Subcommand};

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Trim poly-A tails from FLNC reads. (tama_flnc_polya_cleanup)
    Polya {
        /// FLNC FASTA file. (`-f`)
        #[arg(short = 'f', long)]
        fasta: std::path::PathBuf,
        /// Output prefix. (`-p`)
        #[arg(short = 'p', long)]
        prefix: String,
        /// Minimum read length to keep. (`-m`)
        #[arg(short = 'm', long, default_value_t = 200)]
        min_length: usize,
    },
}

pub fn run(args: Args) -> anyhow::Result<()> {
    match args.cmd {
        Cmd::Polya {
            fasta,
            prefix,
            min_length,
        } => polya(&fasta, &prefix, min_length),
    }
}

const A_PERCENT_THRESHOLD: f64 = 0.7;
const LENGTH_BIN_SIZE: i64 = 200;

struct Block {
    block_type: &'static str,
    #[allow(dead_code)]
    block_count: i64,
    block_start: i64,
    block_end: i64,
}

/// Load FASTA records preserving order: (description, uppercased sequence).
fn load_fasta_records(path: &std::path::Path) -> anyhow::Result<Vec<(String, String)>> {
    let reader = tama_io::open_reader(path)?;
    let mut out: Vec<(String, String)> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            out.push((h.to_string(), String::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.push_str(line.trim_end());
        }
    }
    for (_, seq) in out.iter_mut() {
        seq.make_ascii_uppercase();
    }
    Ok(out)
}

/// Trim poly-A tails from FLNC reads. Faithful port of `tama_flnc_polya_cleanup`.
fn polya(fasta: &std::path::Path, prefix: &str, min_length: usize) -> anyhow::Result<()> {
    let records = load_fasta_records(fasta)?;

    let mut out = tama_io::create_writer(format!("{prefix}.fa"))?;
    let mut out_tails = tama_io::create_writer(format!("{prefix}_tails.fa"))?;
    let mut out_report = tama_io::create_writer(format!("{prefix}_polya_flnc_report.txt"))?;
    let mut out_filtered = tama_io::create_writer(format!("{prefix}_discarded_reads.txt"))?;
    let mut out_summary = tama_io::create_writer(format!("{prefix}_summary.txt"))?;

    let mut polya_dict: BTreeMap<i64, i64> = BTreeMap::new();
    let mut read_length_dict: BTreeMap<i64, i64> = BTreeMap::new();
    let (mut total, mut pass, mut discarded) = (0i64, 0i64, 0i64);
    let mut longest_trimmed_len = 0i64;
    let mut longest_trimmed_id = String::new();
    let mut longest_polya_length = 0i64;

    for (seq_name, seq_string) in &records {
        let seq: Vec<u8> = seq_string.bytes().collect();
        let len = seq.len() as i64;
        total += 1;

        // build A / non-A blocks walking from the 3' end
        let mut blocks: Vec<Block> = Vec::new();
        let mut block_index = 0i64;
        let mut this_block_flag = "NA";
        let mut block_start = -1i64;
        let mut block_end = -1i64;
        let mut this_block_count = 0i64;
        for i in 0..len {
            block_end += 1;
            this_block_count += 1;
            let last_index = (len - 1 - i) as usize;
            let is_a = seq[last_index] == b'A';
            if this_block_flag == "NA" {
                this_block_flag = if is_a { "a_block" } else { "not_a_block" };
                block_start = 0;
                this_block_count = 0;
                continue;
            }
            if is_a && this_block_flag == "a_block" {
                continue;
            }
            if !is_a && this_block_flag == "not_a_block" {
                continue;
            }
            // transition
            blocks.push(Block {
                block_type: this_block_flag,
                block_count: this_block_count,
                block_start,
                block_end,
            });
            block_index += 1;
            this_block_count = 0;
            block_start = block_end;
            this_block_flag = if is_a { "a_block" } else { "not_a_block" };
        }

        // find the tail index
        let mut tail_index = 0i64;
        let mut this_block_index = -1i64;
        let mut stop = false;
        let mut polya_count = 0i64;
        let mut tail_count = 0i64;
        let mut this_a_percent = 1.0f64;
        let mut simple_a_tail_index = 0i64;
        let mut block_overrun = false;

        while !stop {
            this_block_index += 1;
            let b = match get_block(&blocks, this_block_index) {
                Some(b) => b,
                None => {
                    stop = true;
                    block_overrun = true;
                    continue;
                }
            };
            let (btype, bcount, bstart, bend) =
                (b.block_type, b.block_count, b.block_start, b.block_end);
            if this_block_index > 0 {
                this_a_percent = polya_count as f64 / tail_count as f64;
            }
            tail_count += bcount;

            if btype == "a_block" {
                if this_block_index == 0 {
                    simple_a_tail_index = -bend;
                }
                match get_block(&blocks, this_block_index + 1) {
                    None => {
                        tail_index = -bend;
                        stop = true;
                        continue;
                    }
                    Some(next) => {
                        polya_count += bcount;
                        if next.block_count > 2 {
                            tail_index = -bend;
                            stop = true;
                        }
                    }
                }
            }
            if btype == "not_a_block" {
                match get_block(&blocks, this_block_index + 1) {
                    None => {
                        tail_index = -bstart;
                        stop = true;
                        continue;
                    }
                    Some(next) => {
                        if bcount > 2 {
                            tail_index = -bstart;
                            stop = true;
                            continue;
                        } else if next.block_count < 2 {
                            if get_block(&blocks, this_block_index + 3).is_none() {
                                tail_index = -bstart;
                                stop = true;
                                continue;
                            }
                            let b_next = get_block(&blocks, this_block_index + 2).unwrap();
                            let c_next = get_block(&blocks, this_block_index + 3).unwrap();
                            if b_next.block_count > 1 || c_next.block_count < 3 {
                                tail_index = -bstart;
                                stop = true;
                            }
                        }
                    }
                }
            }

            if this_block_index > 2 && this_a_percent < A_PERCENT_THRESHOLD {
                tail_index = if btype == "a_block" { -bend } else { -bstart };
                stop = true;
            }
            if this_block_index >= block_index {
                stop = true;
            }
        }

        // slice the tail and verify A content
        let tail_start = clamp_neg_index(tail_index, len);
        let tail_string = &seq[tail_start..];
        let a_count = tail_string.iter().filter(|&&c| c == b'A').count();
        let a_percent = if tail_string.is_empty() {
            0.0
        } else {
            a_count as f64 / tail_string.len() as f64
        };
        if a_percent < A_PERCENT_THRESHOLD {
            tail_index = if simple_a_tail_index != 0 {
                simple_a_tail_index
            } else {
                0
            };
        }

        let (trim, tail): (&[u8], &[u8]) = if tail_index > -1 {
            (&seq[..], &[])
        } else {
            let cut = clamp_neg_index(tail_index, len);
            (&seq[..cut], &seq[cut..])
        };
        let trim_string = String::from_utf8_lossy(trim);
        let tail_string_out = String::from_utf8_lossy(tail);

        *polya_dict.entry(polya_count).or_insert(0) += 1;

        if trim.len() < min_length {
            writeln!(
                out_filtered,
                ">{seq_name}\tshort\ttrimlen:{}\tprelen:{}",
                trim.len(),
                seq.len()
            )?;
            writeln!(out_filtered, "{seq_string}")?;
            discarded += 1;
        } else if block_overrun {
            writeln!(
                out_filtered,
                ">{seq_name}\toverrun\ttrimlen:{}\tprelen:{}",
                trim.len(),
                seq.len()
            )?;
            writeln!(out_filtered, "{seq_string}")?;
            discarded += 1;
        } else {
            writeln!(out, ">{seq_name}")?;
            writeln!(out, "{trim_string}")?;
            pass += 1;
            let tlen = trim.len() as i64;
            if tlen > longest_trimmed_len {
                longest_trimmed_len = tlen;
                longest_trimmed_id = seq_name.clone();
            }
            let bin_length = (tlen / LENGTH_BIN_SIZE) * LENGTH_BIN_SIZE;
            *read_length_dict.entry(bin_length).or_insert(0) += 1;
        }

        writeln!(out_tails, ">tail_{seq_name}")?;
        writeln!(out_tails, "{tail_string_out}")?;
    }

    // poly-A report + summary counts
    writeln!(out_report, "polya_num\tpolya_num_count")?;
    let (mut zero, mut upto5, mut six_more) = (0i64, 0i64, 0i64);
    for (&polya_num, &count) in &polya_dict {
        writeln!(out_report, "{polya_num}\t{count}")?;
        if polya_num == 0 {
            zero = count;
        } else if polya_num > 0 && polya_num < 6 {
            upto5 += count;
        } else if polya_num > 6 {
            six_more += count;
        }
        if polya_num > longest_polya_length {
            longest_polya_length = polya_num;
        }
    }

    let (mut peak_len, mut peak_count) = (0i64, 0i64);
    for (&bin_length, &count) in &read_length_dict {
        if count > peak_count {
            peak_count = count;
            peak_len = bin_length;
        }
    }

    writeln!(out_summary, "total_read_count:\t{total}")?;
    writeln!(out_summary, "pass_read_count:\t{pass}")?;
    writeln!(out_summary, "discarded_read_count:\t{discarded}")?;
    writeln!(out_summary, "polya_zero_read_count:\t{zero}")?;
    writeln!(out_summary, "polya_upto_five_read_count:\t{upto5}")?;
    writeln!(out_summary, "polya_sixormore_read_count:\t{six_more}")?;
    writeln!(out_summary, "longest_trimmed_read_length:\t{longest_trimmed_len}\tlongest_trimmed_read_id:\t{longest_trimmed_id}")?;
    writeln!(
        out_summary,
        "peak_read_length:\t{peak_len}\tpeak_read_count:\t{peak_count}"
    )?;
    writeln!(out_summary, "longest_polya_length:\t{longest_polya_length}")?;
    Ok(())
}

/// Block at `idx` if present (the original used a dict keyed by block index).
fn get_block(v: &[Block], idx: i64) -> Option<&Block> {
    if idx >= 0 && (idx as usize) < v.len() {
        Some(&v[idx as usize])
    } else {
        None
    }
}

/// Convert a Python negative slice index into a start offset (`seq[idx:]`).
fn clamp_neg_index(idx: i64, len: i64) -> usize {
    if idx >= 0 {
        (idx.min(len)) as usize
    } else {
        (len + idx).max(0) as usize
    }
}
