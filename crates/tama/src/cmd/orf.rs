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
    AddCds {
        /// Blastp parse file. (`-p`)
        #[arg(short = 'p', long)]
        parse: std::path::PathBuf,
        /// Annotation BED file. (`-a`)
        #[arg(short = 'a', long)]
        bed: std::path::PathBuf,
        /// Transcript FASTA. (`-f`)
        #[arg(short = 'f', long)]
        fasta: std::path::PathBuf,
        /// Output BED file. (`-o`)
        #[arg(short = 'o', long)]
        output: std::path::PathBuf,
        /// Include the stop codon in the CDS (`include_stop`). (`-s`)
        #[arg(short = 's', long, default_value = "no_stop_codon")]
        stop: String,
        /// Distance from last SJ to call NMD. (`-d`)
        #[arg(short = 'd', long, default_value_t = 50)]
        sj_dist: i64,
    },
    /// Parse blastp output for ORF selection. (tama_orf_blastp_parser)
    BlastpParse,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    match args.cmd {
        Cmd::Seek { fasta, output } => seek(&fasta, &output),
        Cmd::ExtractCds => Err(super::not_implemented("orf extract-cds")),
        Cmd::AddCds { parse, bed, fasta, output, stop, sj_dist } => {
            add_cds(&parse, &bed, &fasta, &output, &stop, sj_dist)
        }
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

/// Assign CDS/UTR regions and NMD flags to a bed from a blastp parse file.
/// Ports `tama_cds_regions_bed_add`.
fn add_cds(
    parse: &std::path::Path,
    bed: &std::path::Path,
    fasta: &std::path::Path,
    output: &std::path::Path,
    stop: &str,
    sj_dist: i64,
) -> anyhow::Result<()> {
    // blastp parse: trans_id -> [id, frame, nuc_start, nuc_end, prot_start, prot_end, prot_id, match_flag]
    let mut trans_dict: indexmap::IndexMap<String, Vec<String>> = indexmap::IndexMap::new();
    for line in read_lines(parse)? {
        let cols: Vec<String> = line.split('\t').map(String::from).collect();
        trans_dict.insert(cols[0].clone(), cols);
    }
    // transcript lengths from fasta (id = header up to ':' or '(')
    let mut trans_len: indexmap::IndexMap<String, i64> = indexmap::IndexMap::new();
    for (id, seq) in load_fasta(fasta)? {
        let tid = id.split(':').next().unwrap_or(&id);
        let tid = tid.split('(').next().unwrap_or(tid);
        trans_len.insert(tid.to_string(), seq.len() as i64);
    }

    let mut out = tama_io::create_writer(output)?;
    for line in read_lines(bed)? {
        let mut cols: Vec<String> = line.split('\t').map(String::from).collect();
        let trans_id = cols[3].clone();
        let block_list: Vec<i64> = cols[10].split(',').filter(|s| !s.is_empty()).map(|s| s.parse().unwrap()).collect();

        let td = match trans_dict.get(&trans_id) {
            None => {
                cols[6] = "0".into();
                cols[7] = "0".into();
                cols[3] = format!("{};none;missing;no_orf;missing;na", cols[3]);
                writeln!(out, "{}", cols.join("\t"))?;
                continue;
            }
            Some(td) => td.clone(),
        };
        let frame = &td[1];
        let prot_id = &td[6];
        let match_flag = &td[7];
        if match_flag == "missing_nucleotides" {
            cols[6] = "0".into();
            cols[7] = "0".into();
            cols[3] = format!("{};{prot_id};missing;{match_flag};missing;missing", cols[3]);
            writeln!(out, "{}", cols.join("\t"))?;
            continue;
        }

        let trans_start: i64 = cols[1].parse()?;
        let block_start_list: Vec<i64> = cols[11].split(',').filter(|s| !s.is_empty()).map(|s| s.parse().unwrap()).collect();
        let strand = cols[5].chars().next().unwrap_or('+');
        let prot_start: i64 = td[4].parse::<i64>()? - 1;
        let nuc_start: i64 = td[2].parse()?;
        let nuc_end: i64 = td[3].parse()?;
        let tlen = trans_len[&trans_id];

        let (cds_rel_start, cds_rel_end) = if strand == '+' {
            let end = if stop == "include_stop" { nuc_end + 1 } else { nuc_end - 2 };
            (nuc_start, end)
        } else {
            let start = if stop == "include_stop" {
                tlen - (nuc_end + 1)
            } else {
                tlen - (nuc_end - 2)
            };
            (start, tlen - nuc_start)
        };

        let mut block_sum = 0i64;
        let (mut exon_cds_start, mut exon_cds_end) = (0usize, 0usize);
        let (mut cds_coord_start, mut cds_coord_end) = (0i64, 0i64);
        let mut exon_start_list = Vec::new();
        let mut exon_end_list = Vec::new();
        for (i, &block_size) in block_list.iter().enumerate() {
            let prev = block_sum;
            block_sum += block_size;
            let es = trans_start + block_start_list[i];
            exon_start_list.push(es);
            exon_end_list.push(es + block_size);
            if cds_rel_start >= prev && cds_rel_start < block_sum {
                exon_cds_start = i;
                cds_coord_start = trans_start + block_start_list[i] + cds_rel_start - prev;
            }
            if cds_rel_end >= prev && cds_rel_end <= block_sum {
                exon_cds_end = i;
                cds_coord_end = trans_start + block_start_list[i] + cds_rel_end - prev;
            }
        }

        let exon_nums = block_list.len();
        let mut nmd_flag = "prot_ok".to_string();
        if strand == '+' {
            if exon_nums > exon_cds_end + 1 {
                let stop_dist = exon_start_list[exon_nums - 1] - cds_coord_end;
                if stop_dist > sj_dist {
                    nmd_flag = format!("NMD{}", exon_nums - (exon_cds_end + 1));
                }
            }
        } else if exon_cds_start > 0 {
            let stop_dist = cds_coord_start - exon_end_list[0];
            if stop_dist > sj_dist {
                nmd_flag = format!("NMD{exon_cds_start}");
            }
        }

        cols[6] = cds_coord_start.to_string();
        cols[7] = cds_coord_end.to_string();
        let degrade = if prot_start == 0 { "5prime_degrade" } else { "full_length" };
        cols[3] = format!("{};{prot_id};{degrade};{match_flag};{nmd_flag};{frame}", cols[3]);
        writeln!(out, "{}", cols.join("\t"))?;
    }
    Ok(())
}

fn read_lines(path: &std::path::Path) -> anyhow::Result<Vec<String>> {
    let reader = tama_io::open_reader(path)?;
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if !line.trim().is_empty() {
            out.push(line);
        }
    }
    Ok(out)
}

/// Load a FASTA as (id, uppercased seq) in order; id = first whitespace token.
fn load_fasta(path: &std::path::Path) -> anyhow::Result<Vec<(String, String)>> {
    let reader = tama_io::open_reader(path)?;
    let mut out: Vec<(String, String)> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            out.push((h.split_whitespace().next().unwrap_or("").to_string(), String::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.push_str(line.trim_end());
        }
    }
    for (_, s) in out.iter_mut() {
        s.make_ascii_uppercase();
    }
    Ok(out)
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
