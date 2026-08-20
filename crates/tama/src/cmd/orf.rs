//! `tama orf` — ORF / NMD prediction tools.

use std::io::{BufRead, Write};

use clap::{Args as ClapArgs, Subcommand};

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Find ORFs in transcript sequences. (tama_orf_seeker)
    Seek {
        /// Transcript FASTA file. (`-f`)
        #[arg(short = 'f', long)]
        fasta: std::path::PathBuf,
        /// Output protein FASTA. (`-o`)
        #[arg(short = 'o', long)]
        output: std::path::PathBuf,
    },
    /// Extract CDS regions from a BED. (tama_bed_extract_cds)
    ExtractCds,
    /// Add CDS regions to a BED. (tama_cds_regions_bed_add)
    AddCds,
    /// Parse blastp output for ORF selection. (tama_orf_blastp_parser)
    BlastpParse,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    match args.cmd {
        Cmd::Seek { fasta, output } => seek(&fasta, &output),
        Cmd::ExtractCds => Err(super::not_implemented("orf extract-cds")),
        Cmd::AddCds => Err(super::not_implemented("orf add-cds")),
        Cmd::BlastpParse => Err(super::not_implemented("orf blastp-parse")),
    }
}

/// Standard codon table used by tama_orf_seeker (`$` = stop).
fn codon(c: &[u8]) -> Option<char> {
    let s: [u8; 3] = [c[0], c[1], c[2]];
    Some(match &s {
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        b"GAT" | b"GAC" => 'D',
        b"GAA" | b"GAG" => 'E',
        b"GGT" | b"GGC" | b"GGG" | b"GGA" => 'G',
        b"TTT" | b"TTC" => 'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
        b"TAT" | b"TAC" => 'Y',
        b"TAA" | b"TAG" | b"TGA" => '$',
        b"TGT" | b"TGC" => 'C',
        b"TGG" => 'W',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        b"CAT" | b"CAC" => 'H',
        b"CAA" | b"CAG" => 'Q',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
        b"ATT" | b"ATC" | b"ATA" => 'I',
        b"ATG" => 'M',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        b"AAT" | b"AAC" => 'N',
        b"AAA" | b"AAG" => 'K',
        _ => return None,
    })
}

struct Orf {
    seq: String,
    a_start: i64,
    a_end: i64,
    n_start: i64,
    n_end: i64,
    frame: i64,
    length: usize,
    start_codon: char,
}

impl Orf {
    fn new(seq: String, a_start: i64, a_end: i64, n_start: i64, n_end: i64, frame: i64) -> Orf {
        let start_codon = seq.chars().next().unwrap_or(' ');
        let length = seq.chars().count();
        Orf { seq, a_start, a_end, n_start, n_end, frame, length, start_codon }
    }
    fn orf_id(&self) -> String {
        format!(
            "{}_{}_{}_{}_{}_{}_{}",
            self.frame, self.a_start, self.a_end, self.n_start, self.n_end, self.length,
            self.start_codon
        )
    }
}

/// Translate one frame; returns (protein, has_N, nuc positions).
fn frame_iterate(seq: &[u8], frame: usize) -> (Vec<char>, bool, Vec<i64>) {
    let mut prot = Vec::new();
    let mut n_pos = Vec::new();
    let mut pos = frame - 1;
    let n = seq.len();
    if n < 2 {
        return (prot, false, n_pos);
    }
    while pos < n - 2 {
        n_pos.push(pos as i64);
        match codon(&seq[pos..pos + 3]) {
            Some(a) => prot.push(a),
            None => return (prot, true, n_pos),
        }
        pos += 3;
    }
    (prot, false, n_pos)
}

/// Find ORFs in a translated frame. Faithful port of `orf_seeker`.
fn orf_seeker(prot: &[char], frame: i64, n_pos: &[i64]) -> Vec<Orf> {
    let mut orfs: std::collections::BTreeMap<String, Orf> = Default::default();

    // leading (possibly 5'-degraded) ORF: from position 0 to the first stop
    let mut orf_seq = String::new();
    for (i, &aa) in prot.iter().enumerate() {
        if aa == '$' {
            if i > 0 {
                let orf = Orf::new(orf_seq.clone(), 1, i as i64, n_pos[0], n_pos[i] + 2, frame);
                orfs.insert(orf.orf_id(), orf);
            }
            break;
        }
        orf_seq.push(aa);
    }

    // ORFs starting with M
    orf_seq.clear();
    let mut orf_start = 0usize;
    for (i, &aa) in prot.iter().enumerate() {
        if orf_seq.is_empty() {
            if aa == 'M' {
                orf_start = i;
                orf_seq.push(aa);
            }
            continue;
        }
        if aa == '$' {
            if !orf_seq.is_empty() {
                let orf = Orf::new(
                    orf_seq.clone(),
                    orf_start as i64 + 1,
                    i as i64,
                    n_pos[orf_start],
                    n_pos[i] + 2,
                    frame,
                );
                orfs.insert(orf.orf_id(), orf);
                orf_seq.clear();
                orf_start = i + 1;
            }
            continue;
        }
        orf_seq.push(aa);
    }

    orfs.into_values().collect()
}

/// Take the ORFs at the four largest lengths across all frames. Port of
/// `sort_orf_list` (descending length; ties keep frame/orf_id order).
fn top_orfs(mut all: Vec<Orf>) -> Vec<Orf> {
    use indexmap::IndexMap;
    let mut by_len: IndexMap<usize, Vec<Orf>> = IndexMap::new();
    for orf in all.drain(..) {
        by_len.entry(orf.length).or_default().push(orf);
    }
    let mut lengths: Vec<usize> = by_len.keys().copied().collect();
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let mut out = Vec::new();
    for &len in lengths.iter().take(4) {
        out.extend(by_len.shift_remove(&len).unwrap());
    }
    out
}

/// Find ORFs in transcript sequences. Ports `tama_orf_seeker`.
fn seek(fasta: &std::path::Path, output: &std::path::Path) -> anyhow::Result<()> {
    let reader = tama_io::open_reader(fasta)?;
    let mut out = tama_io::create_writer(output)?;

    // parse FASTA (id = first token of header)
    let mut records: Vec<(String, String)> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            let id = h.split_whitespace().next().unwrap_or("").to_string();
            records.push((id, String::new()));
        } else if let Some(last) = records.last_mut() {
            last.1.push_str(line.trim_end());
        }
    }

    for (trans_id, seq_str) in &records {
        let seq: Vec<u8> = seq_str.to_uppercase().into_bytes();
        let frames: Vec<(Vec<char>, bool, Vec<i64>)> =
            (1..=3).map(|f| frame_iterate(&seq, f)).collect();
        if frames.iter().any(|(_, n_flag, _)| *n_flag) {
            writeln!(out, ">{trans_id}:missing_nucleotides")?;
            continue;
        }
        let mut all = Vec::new();
        for (f, (prot, _, n_pos)) in frames.iter().enumerate() {
            all.extend(orf_seeker(prot, f as i64 + 1, n_pos));
        }
        for orf in top_orfs(all) {
            writeln!(
                out,
                ">{trans_id}:F{}:{}:{}:{}:{}:{}:{}",
                orf.frame, orf.n_start, orf.n_end, orf.a_start, orf.a_end, orf.length,
                orf.start_codon
            )?;
            writeln!(out, "{}", orf.seq)?;
        }
    }
    Ok(())
}
