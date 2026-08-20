//! `tama format` — format conversion tools.

use clap::{Args as ClapArgs, Subcommand, ValueEnum};

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
    Bed2gtf,
    /// Convert TAMA BED to GTF with ORF/NMD CDS. (tama_convert_bed_gtf_ensembl_orf_nmd)
    Bed2gtfOrf,
    /// Convert Nanopore FASTQ to FASTA. (tama_convert_nanopore_fastq_fasta)
    Fastq2fasta,
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
    /// Filter transcript models by ID list. (tama_format_id_filter)
    IdFilter,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let name = match args.cmd {
        Cmd::Bed2gtf => "format bed2gtf",
        Cmd::Bed2gtfOrf => "format bed2gtf-orf",
        Cmd::Fastq2fasta => "format fastq2fasta",
        Cmd::Gtf2bed { .. } => "format gtf2bed",
        Cmd::Gff2bed { .. } => "format gff2bed",
        Cmd::IdFilter => "format id-filter",
    };
    Err(super::not_implemented(name))
}
