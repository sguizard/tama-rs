//! `tama format` — format conversion tools.

use std::io::Write;

use anyhow::bail;
use clap::{Args as ClapArgs, Subcommand, ValueEnum};
use indexmap::IndexMap;

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Clone, ValueEnum)]
pub enum GtfSource {
    Ensembl,
    Ncbi,
    Stringtie,
}

#[derive(Clone, ValueEnum)]
pub enum GffSource {
    Cupcake,
    Liftoff,
}

#[derive(Subcommand)]
enum Cmd {
    /// Convert TAMA BED to Ensembl-style GTF (no CDS). (tama_convert_bed_gtf_ensembl_no_cds)
    Bed2gtf {
        /// Input TAMA BED file.
        bed: std::path::PathBuf,
        /// Output GTF file.
        output: std::path::PathBuf,
    },
    /// Convert TAMA BED to GTF with ORF/NMD CDS. (tama_convert_bed_gtf_ensembl_orf_nmd)
    Bed2gtfOrf,
    /// Convert Nanopore FASTQ to FASTA. (tama_convert_nanopore_fastq_fasta)
    Fastq2fasta {
        /// Input FASTQ file.
        fastq: std::path::PathBuf,
        /// Output FASTA file.
        output: std::path::PathBuf,
    },
    /// Convert a GTF to TAMA BED12. (tama_format_gtf_to_bed12_*)
    Gtf2bed {
        #[arg(long)]
        source: GtfSource,
    },
    /// Convert a GFF to TAMA BED12. (tama_format_gff_to_bed12_*)
    Gff2bed {
        #[arg(long)]
        source: GffSource,
    },
    /// Restructure/filter BED ID fields. (tama_format_id_filter)
    IdFilter {
        /// BED file. (`-b`)
        #[arg(short = 'b', long)]
        bed: std::path::PathBuf,
        /// Output file. (`-o`)
        #[arg(short = 'o', long)]
        output: std::path::PathBuf,
        /// Filter level: `none` or `only_match`. (`-f`)
        #[arg(short = 'f', long, default_value = "none")]
        filter: String,
        /// Sub-field method: `ensembl_merge`, `ensembl_orf`, or `custom`. (`-s`)
        #[arg(short = 's', long, default_value = "ensembl_merge")]
        method: String,
        /// Custom reshuffle parameter, e.g. `3,4,1,2`. (`-r`)
        #[arg(short = 'r', long, default_value = "none")]
        reshuffle: String,
        /// Sub-field delimiters. (`-d`)
        #[arg(short = 'd', long, default_value = ";")]
        delim: String,
    },
}

pub fn run(args: Args) -> anyhow::Result<()> {
    match args.cmd {
        Cmd::Bed2gtf { bed, output } => bed2gtf(&bed, &output),
        Cmd::Fastq2fasta { fastq, output } => fastq2fasta(&fastq, &output),
        Cmd::IdFilter {
            bed,
            output,
            filter,
            method,
            reshuffle,
            delim,
        } => id_filter(&bed, &output, &filter, &method, &reshuffle, &delim),
        Cmd::Bed2gtfOrf => Err(super::not_implemented("format bed2gtf-orf")),
        Cmd::Gtf2bed { .. } => Err(super::not_implemented("format gtf2bed")),
        Cmd::Gff2bed { .. } => Err(super::not_implemented("format gff2bed")),
    }
}

const GTF_SOURCE: &str = "PBRI";

/// Convert a TAMA BED to Ensembl-style GTF. Ports `tama_convert_bed_gtf_ensembl_no_cds`.
fn bed2gtf(bed: &std::path::Path, output: &std::path::Path) -> anyhow::Result<()> {
    let reader = tama_io::open_reader(bed)?;
    let transcripts = tama_core::bed::read_bed(reader)?;

    // preserve gene order (first appearance) and transcript order (bed order)
    let mut gene_order: Vec<String> = Vec::new();
    let mut gene_trans: IndexMap<String, Vec<usize>> = IndexMap::new();
    for (i, t) in transcripts.iter().enumerate() {
        let gene = t.gene_id().to_string();
        if !gene_trans.contains_key(&gene) {
            gene_order.push(gene.clone());
        }
        gene_trans.entry(gene).or_default().push(i);
    }

    let mut out = tama_io::create_writer(output)?;
    for gene in &gene_order {
        let idxs = &gene_trans[gene];
        let first = &transcripts[idxs[0]];
        let chrom = &first.scaffold;
        let strand = first.strand;
        let min_start = idxs.iter().map(|&i| transcripts[i].trans_start).min().unwrap();
        let max_end = idxs.iter().map(|&i| transcripts[i].trans_end).max().unwrap();

        writeln!(
            out,
            "{chrom}\t{GTF_SOURCE}\tgene\t{}\t{max_end}\t.\t{strand}\t.\tgene_id \"{gene}\";",
            min_start + 1
        )?;

        for &i in idxs {
            let t = &transcripts[i];
            let trans_id = t.trans_id();
            writeln!(
                out,
                "{chrom}\t{GTF_SOURCE}\ttranscript\t{}\t{}\t.\t{strand}\t.\tgene_id \"{gene}\"; transcript_id \"{trans_id}\"; uniq_trans_id \"{trans_id}\";",
                t.trans_start + 1,
                t.trans_end
            )?;
            let n = t.num_exons();
            for k in 0..n {
                let e_num = if strand == '+' { k + 1 } else { n - k };
                writeln!(
                    out,
                    "{chrom}\t{GTF_SOURCE}\texon\t{}\t{}\t.\t{strand}\t.\tgene_id \"{gene}\"; transcript_id \"{trans_id}\"; exon_number \"{e_num}\"; uniq_trans_id \"{trans_id}\";",
                    t.exon_start_list[k] + 1,
                    t.exon_end_list[k]
                )?;
            }
        }
    }
    Ok(())
}

/// Convert a 4-line-per-record FASTQ to FASTA. Ports `tama_convert_nanopore_fastq_fasta`.
fn fastq2fasta(fastq: &std::path::Path, output: &std::path::Path) -> anyhow::Result<()> {
    use std::io::BufRead;
    let reader = tama_io::open_reader(fastq)?;
    let mut out = tama_io::create_writer(output)?;
    let mut count = 0;
    let mut header = String::new();
    for line in reader.lines() {
        let line = line?;
        count += 1;
        if count == 1 && line.starts_with('@') {
            header = format!(">{}", &line[1..]);
        } else if count == 2 {
            writeln!(out, "{header}")?;
            writeln!(out, "{line}")?;
        }
        if count == 4 {
            count = 0;
        }
    }
    Ok(())
}

/// Restructure/filter BED ID fields. Ports `tama_format_id_filter`.
fn id_filter(
    bed: &std::path::Path,
    output: &std::path::Path,
    filter: &str,
    method: &str,
    reshuffle: &str,
    delim: &str,
) -> anyhow::Result<()> {
    use std::io::BufRead;
    let reader = tama_io::open_reader(bed)?;
    let mut out = tama_io::create_writer(output)?;
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let mut cols: Vec<String> = line.split('\t').map(String::from).collect();
        if cols.len() < 12 {
            continue;
        }
        // `None` = filtered out (only_match with no match).
        if let Some(new_id) = id_parse(&cols[3], filter, method, reshuffle, delim)? {
            cols[3] = new_id;
            writeln!(out, "{}", cols.join("\t"))?;
        }
    }
    Ok(())
}

/// Returns `Some(new_id_line)` to keep the record, or `None` to drop it.
fn id_parse(
    id_line: &str,
    filter: &str,
    method: &str,
    reshuffle: &str,
    delim: &str,
) -> anyhow::Result<Option<String>> {
    match method {
        "ensembl_merge" => {
            let parts: Vec<&str> = id_line.split(delim).collect();
            let tama_gene = parts[0];
            let tama_trans = parts.get(1).copied().unwrap_or("");
            let mut ens_gene = String::new();
            let mut ens_trans = String::new();
            let mut other: Vec<&str> = Vec::new();
            if parts.len() > 3 {
                ens_gene = parts[2].to_string();
                ens_trans = parts[3].to_string();
                if !ens_gene.contains("ENS") {
                    bail!("Ensembl gene ID is not right: {ens_gene}");
                }
                if !ens_trans.contains("ENS") {
                    bail!("Ensembl transcript ID is not right: {ens_trans}");
                }
                other = parts[4..].to_vec();
            }
            match filter {
                "only_match" => {
                    if !ens_gene.contains("ENS") {
                        return Ok(None);
                    }
                }
                "none" => {
                    if ens_gene.is_empty() {
                        ens_gene = tama_gene.to_string();
                    }
                    if ens_trans.is_empty() {
                        ens_trans = tama_trans.to_string();
                    }
                }
                other => bail!("invalid filter level {other:?}"),
            }
            let mut new = vec![ens_gene, ens_trans, tama_gene.to_string(), tama_trans.to_string()];
            new.extend(other.iter().map(|s| s.to_string()));
            Ok(Some(new.join(";")))
        }
        "ensembl_orf" => {
            let parts: Vec<&str> = id_line.split(delim).collect();
            let tama_gene = parts[0];
            let tama_trans = parts.get(1).copied().unwrap_or("");
            let mut ens_gene = String::new();
            let mut ens_trans = String::new();
            let mut other: Vec<&str> = Vec::new();
            if parts.len() > 3 {
                let ens_id_split: Vec<&str> = parts[2].split(',').collect();
                if ens_id_split.len() > 1 {
                    ens_gene = ens_id_split[0].to_string();
                    ens_trans = ens_id_split[1].to_string();
                    if !ens_gene.contains("ENS") {
                        bail!("Ensembl gene ID is not right: {ens_gene}");
                    }
                    if !ens_trans.contains("ENS") {
                        bail!("Ensembl transcript ID is not right: {ens_trans}");
                    }
                } else {
                    ens_gene = tama_gene.to_string();
                    ens_trans = tama_trans.to_string();
                }
                other = parts[3..].to_vec();
            }
            if filter == "only_match" && !ens_gene.contains("ENS") {
                return Ok(None);
            }
            let mut new = vec![ens_gene, ens_trans, tama_gene.to_string(), tama_trans.to_string()];
            new.extend(other.iter().map(|s| s.to_string()));
            Ok(Some(new.join(";")))
        }
        "custom" => {
            // split id_line on any of the delimiter characters
            let delims: Vec<char> = delim.chars().collect();
            let fields: Vec<&str> = id_line.split(|c| delims.contains(&c)).collect();
            let mut new = Vec::new();
            for idx in reshuffle.split(',') {
                let i: usize = idx.trim().parse::<usize>()?;
                new.push(fields.get(i - 1).copied().unwrap_or("").to_string());
            }
            Ok(Some(new.join(";")))
        }
        other => bail!("invalid subfield method {other:?}"),
    }
}
