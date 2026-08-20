//! `tama split` — file splitting tools.

use std::io::{BufRead, Write};

use clap::{Args as ClapArgs, Subcommand};
use indexmap::IndexMap;

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Split a FASTA into chunks. (tama_fasta_splitter)
    Fasta {
        /// Input FASTA file.
        fasta: std::path::PathBuf,
        /// Output prefix (files are `<prefix>_<n>.fa`).
        prefix: String,
        /// Number of chunks to split into.
        split_number: usize,
    },
    /// Split a mapped SAM by chromosome. (tama_mapped_sam_splitter)
    Sam {
        /// Input SAM file.
        sam: std::path::PathBuf,
        /// Number of output files.
        num_files: usize,
        /// Output prefix (files are `<prefix>_<n>.sam`).
        prefix: String,
    },
}

pub fn run(args: Args) -> anyhow::Result<()> {
    match args.cmd {
        Cmd::Fasta { fasta, prefix, split_number } => split_fasta(&fasta, &prefix, split_number),
        Cmd::Sam { sam, num_files, prefix } => split_sam(&sam, num_files, &prefix),
    }
}

/// Split a FASTA into `split_number` roughly-equal files by sequence count.
/// Ports `tama_fasta_splitter`.
fn split_fasta(fasta: &std::path::Path, prefix: &str, split_number: usize) -> anyhow::Result<()> {
    let reader = tama_io::open_reader(fasta)?;
    // preserve order; each header maps to its full record (header + seq lines)
    let mut records: Vec<Vec<String>> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            records.push(vec![line]);
        } else if let Some(last) = records.last_mut() {
            last.push(line);
        }
    }

    let seq_count = records.len();
    if split_number == 0 {
        anyhow::bail!("split number must be > 0");
    }
    // ceil(seq_count / split_number)
    let split_size = seq_count.div_ceil(split_number).max(1);

    let mut split_index = 0;
    let mut out: Option<Box<dyn Write>> = None;
    for (i, rec) in records.iter().enumerate() {
        if i % split_size == 0 {
            split_index += 1;
            out = Some(tama_io::create_writer(format!("{prefix}_{split_index}.fa"))?);
        }
        let w = out.as_mut().unwrap();
        for item in rec {
            writeln!(w, "{item}")?;
        }
    }
    Ok(())
}

/// Split a mapped SAM into `num_files` files, keeping whole chromosomes together.
/// Ports `tama_mapped_sam_splitter`.
fn split_sam(sam: &std::path::Path, num_files: usize, prefix: &str) -> anyhow::Result<()> {
    if num_files == 0 {
        anyhow::bail!("number of files must be > 0");
    }
    let reader = tama_io::open_reader(sam)?;
    let mut headers: Vec<String> = Vec::new();
    let mut chrom_order: Vec<String> = Vec::new();
    let mut chrom_reads: IndexMap<String, Vec<String>> = IndexMap::new();
    let mut sam_count = 0usize;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            headers.push(line);
            continue;
        }
        if line.is_empty() {
            continue;
        }
        sam_count += 1;
        let scaff = line.split('\t').nth(2).unwrap_or("").to_string();
        if !chrom_reads.contains_key(&scaff) {
            chrom_order.push(scaff.clone());
        }
        chrom_reads.entry(scaff).or_default().push(line);
    }

    let num_file_reads = sam_count / num_files; // py2 integer division

    let mut file_count = 1;
    let mut file_read_count = 0usize;
    let mut total = 0usize;
    let mut chrom_idx = 0usize;

    let mut out = tama_io::create_writer(format!("{prefix}_{file_count}.sam"))?;
    for h in &headers {
        writeln!(out, "{h}")?;
    }

    while total < sam_count {
        let chrom = &chrom_order[chrom_idx];
        let chrom_numreads = chrom_reads[chrom].len();
        let prev_file_read_count = file_read_count;
        file_read_count += chrom_numreads;
        chrom_idx += 1;

        if file_read_count > num_file_reads && prev_file_read_count > 0 {
            file_read_count = chrom_numreads;
            file_count += 1;
            out = tama_io::create_writer(format!("{prefix}_{file_count}.sam"))?;
            for h in &headers {
                writeln!(out, "{h}")?;
            }
        }

        for read_line in &chrom_reads[chrom] {
            writeln!(out, "{read_line}")?;
            total += 1;
        }
    }
    Ok(())
}
