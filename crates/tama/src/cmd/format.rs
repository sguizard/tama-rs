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
    Bed2gtfOrf {
        /// Input ORF/NMD BED file.
        bed: std::path::PathBuf,
        /// Output GTF file.
        output: std::path::PathBuf,
    },
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
        /// Input GTF file.
        gtf: std::path::PathBuf,
        /// Output BED file.
        output: std::path::PathBuf,
    },
    /// Convert a GFF to TAMA BED12. (tama_format_gff_to_bed12_*)
    Gff2bed {
        #[arg(long)]
        source: GffSource,
        /// Input GFF file.
        gff: std::path::PathBuf,
        /// Output BED file.
        output: std::path::PathBuf,
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
        Cmd::Bed2gtfOrf { bed, output } => bed2gtf_orf(&bed, &output),
        Cmd::Gtf2bed {
            source,
            gtf,
            output,
        } => match source {
            GtfSource::Stringtie => gtf2bed_stringtie(&gtf, &output),
            GtfSource::Ensembl => gtf2bed_ensembl(&gtf, &output),
            GtfSource::Ncbi => gtf2bed_ncbi(&gtf, &output),
        },
        Cmd::Gff2bed {
            source,
            gff,
            output,
        } => match source {
            GffSource::Cupcake => gff2bed_cupcake(&gff, &output),
            GffSource::Liftoff => gff2bed_liftoff(&gff, &output),
        },
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
        let min_start = idxs
            .iter()
            .map(|&i| transcripts[i].trans_start)
            .min()
            .unwrap();
        let max_end = idxs
            .iter()
            .map(|&i| transcripts[i].trans_end)
            .max()
            .unwrap();

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

/// A transcript for ORF-GTF conversion.
struct OrfTx {
    chrom: String,
    t_start: i64, // 1-based
    t_end: i64,
    gene_id: String,
    trans_id: String,
    strand: char,
    num_exons: usize,
    starts: Vec<i64>, // 1-based exon starts
    ends: Vec<i64>,
    cds_start: i64, // 1-based
    cds_end: i64,
    prot_id: String,
    degrade: String,
    match_flag: String,
    nmd: String,
}

/// Convert an ORF/NMD BED to GTF with CDS/UTR/codon features. Ports
/// `tama_convert_bed_gtf_ensembl_orf_nmd`.
fn bed2gtf_orf(bed: &std::path::Path, output: &std::path::Path) -> anyhow::Result<()> {
    use std::io::BufRead;
    const SOURCE: &str = "PBRI";

    let anno = |pairs: &[(&str, &str)]| -> String {
        let mut parts: Vec<String> = pairs.iter().map(|(k, v)| format!("{k} \"{v}\"")).collect();
        parts.push(String::new());
        parts.join("; ")
    };

    let mut gene_order: Vec<String> = Vec::new();
    let mut gene_trans: IndexMap<String, Vec<String>> = IndexMap::new();
    let mut trans: IndexMap<String, OrfTx> = IndexMap::new();

    let reader = tama_io::open_reader(bed)?;
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        let id_split: Vec<&str> = c[3].split(';').collect();
        let (gene_id, trans_id) = (id_split[0].to_string(), id_split[1].to_string());
        let t0: i64 = c[1].parse()?;
        let sizes: Vec<i64> = c[10]
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().unwrap())
            .collect();
        let rels: Vec<i64> = c[11]
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().unwrap())
            .collect();
        let starts: Vec<i64> = rels.iter().map(|r| t0 + r + 1).collect();
        // calc_end: start + block_size - 1 (bed 0-based to gtf 1-based)
        let ends: Vec<i64> = starts.iter().zip(&sizes).map(|(s, z)| s + z - 1).collect();
        let tx = OrfTx {
            chrom: c[0].to_string(),
            t_start: t0 + 1,
            t_end: c[2].parse()?,
            gene_id: gene_id.clone(),
            trans_id: trans_id.clone(),
            strand: c[5].chars().next().unwrap_or('+'),
            num_exons: c[9].parse()?,
            starts,
            ends,
            cds_start: c[6].parse::<i64>()? + 1,
            cds_end: c[7].parse()?,
            prot_id: id_split.get(2).unwrap_or(&"").to_string(),
            degrade: id_split.get(3).unwrap_or(&"").to_string(),
            match_flag: id_split.get(4).unwrap_or(&"").to_string(),
            nmd: id_split.get(5).unwrap_or(&"").to_string(),
        };
        if !gene_trans.contains_key(&gene_id) {
            gene_order.push(gene_id.clone());
        }
        gene_trans
            .entry(gene_id)
            .or_default()
            .push(trans_id.clone());
        trans.insert(trans_id, tx);
    }

    let mut out = tama_io::create_writer(output)?;
    for gene_id in &gene_order {
        let tids = &gene_trans[gene_id];
        let first = &trans[&tids[0]];
        let chrom = first.chrom.clone();
        let strand = first.strand;
        let min_start = tids.iter().map(|t| trans[t].t_start).min().unwrap();
        let max_end = tids.iter().map(|t| trans[t].t_end).max().unwrap();
        writeln!(
            out,
            "{chrom}\t{SOURCE}\tgene\t{min_start}\t{max_end}\t.\t{strand}\t.\t{}",
            anno(&[("gene_id", gene_id), ("gene_source", "tama")])
        )?;

        for tid in tids {
            let tx = &trans[tid];
            write_orf_trans(&mut out, tx, SOURCE, &anno)?;
        }
    }
    Ok(())
}

#[allow(clippy::type_complexity)]
fn write_orf_trans(
    out: &mut Box<dyn Write>,
    tx: &OrfTx,
    source: &str,
    anno: &dyn Fn(&[(&str, &str)]) -> String,
) -> anyhow::Result<()> {
    let g = &tx.gene_id;
    let t = &tx.trans_id;
    writeln!(
        out,
        "{}\t{source}\ttranscript\t{}\t{}\t.\t{}\t.\t{}",
        tx.chrom,
        tx.t_start,
        tx.t_end,
        tx.strand,
        anno(&[
            ("gene_id", g),
            ("transcript_id", t),
            ("gene_source", "tama"),
            ("transcript_source", "tama")
        ])
    )?;

    let mut five_utr: Vec<(i64, i64, usize)> = Vec::new();
    let mut three_utr: Vec<(i64, i64, usize)> = Vec::new();

    let feature = |out: &mut Box<dyn Write>,
                   region: &str,
                   s: i64,
                   e: i64,
                   e_num: usize|
     -> anyhow::Result<()> {
        writeln!(
            out,
            "{}\t{source}\t{region}\t{s}\t{e}\t.\t{}\t.\t{}",
            tx.chrom,
            tx.strand,
            anno(&[
                ("gene_id", g),
                ("transcript_id", t),
                ("exon_number", &e_num.to_string()),
                ("gene_source", "tama"),
                ("transcript_source", "tama"),
                ("prot_id", &tx.prot_id),
                ("degrade_flag", &tx.degrade),
                ("match_flag", &tx.match_flag),
                ("nmd_flag", &tx.nmd),
            ])
        )?;
        Ok(())
    };

    let no_cds = tx.cds_start == 1 && tx.cds_end == 0;
    for i in 0..tx.num_exons {
        let e_num = i + 1;
        let e_index = if tx.strand == '+' {
            i
        } else {
            tx.num_exons - i - 1
        };
        let e_start = tx.starts[e_index];
        let e_end = tx.ends[e_index];
        feature(out, "exon", e_start, e_end, e_num)?;
        if no_cds {
            continue;
        }

        if tx.strand == '+' {
            if e_start > tx.cds_end {
                three_utr.push((e_start, e_end, e_num));
                continue;
            }
            if tx.cds_start > e_end {
                five_utr.push((e_start, e_end, e_num));
                continue;
            }
            if e_start > tx.cds_start && e_end < tx.cds_end {
                feature(out, "CDS", e_start, e_end, e_num)?;
                continue;
            }
            let (mut cds_e_start, mut cds_e_end) = (e_start, e_end);
            let (mut start_flag, mut stop_flag) = (false, false);
            let (mut sc_s, mut sc_e, mut stc_s, mut stc_e) = (0, 0, 0, 0);
            if e_start < tx.cds_end && e_end >= tx.cds_end {
                stc_s = tx.cds_end - 2;
                stc_e = tx.cds_end;
                cds_e_end = stc_s - 1;
                three_utr.push((stc_e + 1, e_end, e_num));
                stop_flag = true;
            }
            if tx.cds_start > e_start && tx.cds_start <= e_end {
                sc_s = tx.cds_start;
                sc_e = tx.cds_start + 2;
                cds_e_start = sc_s;
                five_utr.push((e_start, tx.cds_start - 1, e_num));
                start_flag = true;
            }
            if tx.cds_start == tx.t_start {
                start_flag = false;
            }
            if cds_e_end >= cds_e_start {
                feature(out, "CDS", cds_e_start, cds_e_end, e_num)?;
            }
            if start_flag {
                feature(out, "start_codon", sc_s, sc_e, e_num)?;
            }
            if stop_flag {
                feature(out, "stop_codon", stc_s, stc_e, e_num)?;
            }
        } else {
            if e_start > tx.cds_end {
                five_utr.push((e_start, e_end, e_num));
                continue;
            }
            if tx.cds_start > e_end {
                three_utr.push((e_start, e_end, e_num));
                continue;
            }
            if e_start > tx.cds_start && e_end < tx.cds_end {
                feature(out, "CDS", e_start, e_end, e_num)?;
                continue;
            }
            let (mut cds_e_start, mut cds_e_end) = (e_start, e_end);
            let (mut start_flag, mut stop_flag) = (false, false);
            let (mut sc_s, mut sc_e, mut stc_s, mut stc_e) = (0, 0, 0, 0);
            if e_start < tx.cds_end && e_end >= tx.cds_end {
                cds_e_end = tx.cds_end;
                sc_s = tx.cds_end - 2;
                sc_e = tx.cds_end;
                five_utr.push((sc_e + 1, e_end, e_num));
                start_flag = true;
            }
            if tx.cds_start > e_start && tx.cds_start <= e_end {
                stc_s = tx.cds_start;
                stc_e = tx.cds_start + 2;
                cds_e_start = stc_e + 1;
                three_utr.push((e_start, stc_s, e_num));
                stop_flag = true;
            }
            if tx.cds_end == tx.t_end {
                start_flag = false;
            }
            if cds_e_end >= cds_e_start {
                feature(out, "CDS", cds_e_start, cds_e_end, e_num)?;
            }
            if start_flag {
                feature(out, "start_codon", sc_s, sc_e, e_num)?;
            }
            if stop_flag {
                feature(out, "stop_codon", stc_s, stc_e, e_num)?;
            }
        }
    }

    for &(s, e, en) in &five_utr {
        if e >= s {
            feature(out, "five_prime_utr", s, e, en)?;
        }
    }
    for &(s, e, en) in &three_utr {
        if e >= s {
            feature(out, "three_prime_utr", s, e, en)?;
        }
    }
    Ok(())
}

/// One attribute lookup: value of the first attribute whose key contains `key`.
fn gtf_attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    for field in attrs.split(';') {
        if field.contains(key) && field.contains('"') {
            return field.split('"').nth(1);
        }
    }
    None
}

/// A transcript accumulated from Ensembl/NCBI GTF feature lines.
#[derive(Default)]
struct EnsTx {
    chrom: String,
    strand: char,
    t_start: i64,
    t_end: i64,
    gene_id: String,
    gene_class: String,
    trans_class: String,
    e_starts: Vec<i64>,
    e_ends: Vec<i64>,
    c_starts: Vec<i64>,
    c_ends: Vec<i64>,
}

/// Convert an Ensembl GTF to TAMA BED12. Ports `tama_format_gtf_to_bed12_ensembl`.
/// CDS thick bounds are `min(cds_start)-1 .. max(cds_end)` (or the transcript
/// bounds when non-coding). The original's codon/UTR checks only emit warnings,
/// so they do not affect output and are omitted.
fn gtf2bed_ensembl(gtf: &std::path::Path, output: &std::path::Path) -> anyhow::Result<()> {
    let lines = read_all_lines(gtf)?;
    let mut order: Vec<String> = Vec::new();
    let mut trans: IndexMap<String, EnsTx> = IndexMap::new();

    // pass 1: transcripts
    for line in &lines {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 || cols[2] != "transcript" {
            continue;
        }
        let trans_id = gtf_attr(cols[8], "transcript_id").unwrap_or("").to_string();
        let tx = EnsTx {
            chrom: cols[0].to_string(),
            strand: cols[6].chars().next().unwrap_or('+'),
            t_start: cols[3].parse::<i64>()? - 1,
            t_end: cols[4].parse()?,
            gene_id: gtf_attr(cols[8], "gene_id").unwrap_or("").to_string(),
            gene_class: gtf_attr(cols[8], "gene_biotype")
                .unwrap_or("NA")
                .to_string(),
            trans_class: gtf_attr(cols[8], "transcript_biotype")
                .unwrap_or("NA")
                .to_string(),
            ..Default::default()
        };
        if !trans.contains_key(&trans_id) {
            order.push(trans_id.clone());
        }
        trans.insert(trans_id, tx);
    }

    // pass 2: exons and CDS
    for line in &lines {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        let feature = cols[2];
        if feature != "exon" && feature != "CDS" {
            continue;
        }
        let trans_id = gtf_attr(cols[8], "transcript_id").unwrap_or("");
        let Some(tx) = trans.get_mut(trans_id) else {
            continue;
        };
        let (start, end): (i64, i64) = (cols[3].parse()?, cols[4].parse()?);
        if feature == "exon" {
            tx.e_starts.push(start - 1);
            tx.e_ends.push(end);
        } else {
            tx.c_starts.push(start);
            tx.c_ends.push(end);
        }
    }

    let mut out = tama_io::create_writer(output)?;
    for trans_id in &order {
        let tx = &trans[trans_id];
        let mut e_starts = tx.e_starts.clone();
        let mut e_ends = tx.e_ends.clone();
        e_starts.sort_unstable();
        e_ends.sort_unstable();
        let (mut blocks, mut rel_starts) = (String::new(), String::new());
        for k in 0..e_starts.len() {
            if k > 0 {
                blocks.push(',');
                rel_starts.push(',');
            }
            blocks.push_str(&(e_ends[k] - e_starts[k]).to_string());
            rel_starts.push_str(&(e_starts[k] - tx.t_start).to_string());
        }
        let (cds_start, cds_end) = if tx.c_starts.is_empty() {
            (tx.t_start, tx.t_start)
        } else {
            let mut cs = tx.c_starts.clone();
            let mut ce = tx.c_ends.clone();
            cs.sort_unstable();
            ce.sort_unstable();
            (cs[0] - 1, *ce.last().unwrap())
        };
        writeln!(
            out,
            "{}\t{}\t{}\t{};{};{};{}\t40\t{}\t{cds_start}\t{cds_end}\t255,0,0\t{}\t{blocks}\t{rel_starts}",
            tx.chrom, tx.t_start, tx.t_end, tx.gene_id, trans_id, tx.gene_class, tx.trans_class,
            tx.strand, e_starts.len()
        )?;
    }
    Ok(())
}

/// Convert an NCBI GTF to TAMA BED12. Ports `tama_format_gtf_to_bed12_ncbi`.
/// Unlike Ensembl, transcripts are derived from `exon` lines (bounds = min/max
/// exon), IDs are `gene;trans`, and `unknown_transcript_1` records are skipped.
fn gtf2bed_ncbi(gtf: &std::path::Path, output: &std::path::Path) -> anyhow::Result<()> {
    let lines = read_all_lines(gtf)?;
    let mut order: Vec<String> = Vec::new();
    let mut trans: IndexMap<String, EnsTx> = IndexMap::new();

    for line in &lines {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        let feature = cols[2];
        if feature != "exon" && feature != "CDS" {
            continue;
        }
        let trans_id = gtf_attr(cols[8], "transcript_id").unwrap_or("");
        if trans_id == "unknown_transcript_1" {
            continue;
        }
        let (start, end): (i64, i64) = (cols[3].parse()?, cols[4].parse()?);
        let tx = trans.entry(trans_id.to_string()).or_insert_with(|| {
            order.push(trans_id.to_string());
            EnsTx {
                chrom: cols[0].to_string(),
                strand: cols[6].chars().next().unwrap_or('+'),
                gene_id: gtf_attr(cols[8], "gene_id").unwrap_or("").to_string(),
                ..Default::default()
            }
        });
        if feature == "exon" {
            tx.e_starts.push(start - 1);
            tx.e_ends.push(end);
        } else {
            tx.c_starts.push(start);
            tx.c_ends.push(end);
        }
    }

    let mut out = tama_io::create_writer(output)?;
    for trans_id in &order {
        let tx = &trans[trans_id];
        let mut e_starts = tx.e_starts.clone();
        let mut e_ends = tx.e_ends.clone();
        e_starts.sort_unstable();
        e_ends.sort_unstable();
        let t_start = e_starts[0];
        let t_end = *e_ends.last().unwrap();
        let (mut blocks, mut rel_starts) = (String::new(), String::new());
        for k in 0..e_starts.len() {
            if k > 0 {
                blocks.push(',');
                rel_starts.push(',');
            }
            blocks.push_str(&(e_ends[k] - e_starts[k]).to_string());
            rel_starts.push_str(&(e_starts[k] - t_start).to_string());
        }
        let (cds_start, cds_end) = if tx.c_starts.is_empty() {
            (t_start, t_start)
        } else {
            let mut cs = tx.c_starts.clone();
            let mut ce = tx.c_ends.clone();
            cs.sort_unstable();
            ce.sort_unstable();
            (cs[0] - 1, *ce.last().unwrap())
        };
        writeln!(
            out,
            "{}\t{t_start}\t{t_end}\t{};{}\t40\t{}\t{cds_start}\t{cds_end}\t255,0,0\t{}\t{blocks}\t{rel_starts}",
            tx.chrom, tx.gene_id, trans_id, tx.strand, e_starts.len()
        )?;
    }
    Ok(())
}

fn read_all_lines(path: &std::path::Path) -> anyhow::Result<Vec<String>> {
    use std::io::BufRead;
    let reader = tama_io::open_reader(path)?;
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if !line.starts_with('#') {
            out.push(line);
        }
    }
    Ok(out)
}

/// Convert a StringTie/Cufflinks GTF to TAMA BED12. Ports
/// `tama_format_gtf_to_bed12_stringtie`. Exons are read from lines carrying an
/// `exon_number` attribute and ordered by that number (which must be
/// genomic-ascending, as the original asserts).
fn gtf2bed_stringtie(gtf: &std::path::Path, output: &std::path::Path) -> anyhow::Result<()> {
    use std::io::BufRead;

    #[derive(Default)]
    struct Tx {
        chrom: String,
        strand: char,
        exons: std::collections::BTreeMap<i64, (i64, i64)>, // exon_number -> (start, end)
    }

    let mut gene_order: Vec<String> = Vec::new();
    let mut gene_trans: IndexMap<String, Vec<String>> = IndexMap::new();
    let mut trans: IndexMap<(String, String), Tx> = IndexMap::new();

    let reader = tama_io::open_reader(gtf)?;
    let mut cur_gene = String::new();
    let mut cur_trans = String::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        let chrom = cols[0];
        let e_start: i64 = cols[3].parse()?;
        let e_end: i64 = cols[4].parse()?;
        let strand = cols[6].chars().next().unwrap_or('+');
        for field in cols[8].split(';') {
            let field = field.trim();
            if field.is_empty() || !field.contains('"') {
                continue;
            }
            let id_code = field.split('"').nth(1).unwrap_or("");
            if field.contains("gene_id") {
                cur_gene = id_code.to_string();
                if !gene_trans.contains_key(&cur_gene) {
                    gene_order.push(cur_gene.clone());
                    gene_trans.insert(cur_gene.clone(), Vec::new());
                }
            } else if field.contains("transcript_id") {
                cur_trans = id_code.to_string();
                let key = (cur_gene.clone(), cur_trans.clone());
                if !trans.contains_key(&key) {
                    gene_trans
                        .get_mut(&cur_gene)
                        .unwrap()
                        .push(cur_trans.clone());
                    trans.insert(
                        key,
                        Tx {
                            chrom: chrom.to_string(),
                            strand,
                            ..Default::default()
                        },
                    );
                }
            } else if field.contains("exon_number") {
                let e_num: i64 = id_code.parse()?;
                let tx = trans
                    .get_mut(&(cur_gene.clone(), cur_trans.clone()))
                    .unwrap();
                if tx.exons.insert(e_num, (e_start, e_end)).is_some() {
                    anyhow::bail!("duplicate exon number {e_num}");
                }
            }
        }
    }

    let mut records = Vec::new();
    for gene in &gene_order {
        for trans_id in &gene_trans[gene] {
            let tx = &trans[&(gene.clone(), trans_id.clone())];
            let exons: Vec<(i64, i64)> = tx.exons.values().copied().collect();
            records.push((
                tx.chrom.clone(),
                tx.strand,
                gene.clone(),
                trans_id.clone(),
                exons,
            ));
        }
    }
    write_gtf_bed(&records, output)
}

type GtfBedRecord = (String, char, String, String, Vec<(i64, i64)>);

/// Write records as TAMA BED12. Each record is (chrom, strand, gene, trans,
/// ascending exons). Shared by the GTF/GFF converters; the coordinate math
/// matches the originals (`t_start = first_exon_start - 1`, block = end+1-start).
fn write_gtf_bed(records: &[GtfBedRecord], output: &std::path::Path) -> anyhow::Result<()> {
    let mut out = tama_io::create_writer(output)?;
    for (chrom, strand, gene, trans_id, exons) in records {
        let t_start = exons[0].0 - 1;
        let t_end = exons.last().unwrap().1;
        let (mut blocks, mut rel_starts) = (String::new(), String::new());
        for (k, &(s, e)) in exons.iter().enumerate() {
            if k > 0 {
                blocks.push(',');
                rel_starts.push(',');
            }
            blocks.push_str(&(e + 1 - s).to_string());
            rel_starts.push_str(&(s - t_start - 1).to_string());
        }
        writeln!(
            out,
            "{chrom}\t{t_start}\t{t_end}\t{gene};{trans_id}\t40\t{strand}\t{t_start}\t{t_end}\t255,0,0\t{}\t{blocks}\t{rel_starts}",
            exons.len()
        )?;
    }
    Ok(())
}

/// GFF3 attribute lookup (`key=value;` fields).
fn gff_attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    for field in attrs.split(';') {
        if let Some((k, v)) = field.split_once('=') {
            if k == key {
                return Some(v);
            }
        }
    }
    None
}

/// Convert a Liftoff GFF3 to TAMA BED12. Ports `tama_format_gff_to_bed12_liftoff`.
/// id = `gene;trans;name;feature_type`; CDS thick = `min(cds)-1 .. max(cds)` (or
/// 0..0 if non-coding).
fn gff2bed_liftoff(gff: &std::path::Path, output: &std::path::Path) -> anyhow::Result<()> {
    const TRANS_TYPES: &[&str] = &[
        "transcript",
        "mRNA",
        "ncRNA",
        "tRNA",
        "rRNA",
        "snRNA",
        "miRNA",
        "primary_transcript",
        "snoRNA",
        "V_gene_segment",
        "guide_RNA",
        "C_gene_segment",
        "telomerase_RNA",
        "SRP_RNA",
        "lnc_RNA",
    ];

    #[derive(Default)]
    struct Tx {
        chrom: String,
        strand: char,
        t_start: i64,
        t_end: i64,
        gene_id: String,
        trans_name: String,
        trans_type: String,
        e_starts: Vec<i64>,
        e_ends: Vec<i64>,
        c_starts: Vec<i64>,
        c_ends: Vec<i64>,
    }

    let lines = read_all_lines(gff)?;
    let mut order: Vec<String> = Vec::new();
    let mut trans: IndexMap<String, Tx> = IndexMap::new();

    // pass 1: transcript-like features
    for line in &lines {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 || !TRANS_TYPES.contains(&cols[2]) {
            continue;
        }
        let trans_id = gff_attr(cols[8], "ID").unwrap_or("").to_string();
        let tx = Tx {
            chrom: cols[0].to_string(),
            strand: cols[6].chars().next().unwrap_or('+'),
            t_start: cols[3].parse::<i64>()? - 1,
            t_end: cols[4].parse()?,
            gene_id: gff_attr(cols[8], "Parent").unwrap_or("").to_string(),
            trans_name: gff_attr(cols[8], "Name").unwrap_or("NA").to_string(),
            trans_type: cols[2].to_string(),
            ..Default::default()
        };
        if !trans.contains_key(&trans_id) {
            order.push(trans_id.clone());
        }
        trans.insert(trans_id, tx);
    }

    // pass 2: exon / CDS features
    for line in &lines {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        let feature = cols[2];
        if feature != "exon" && feature != "CDS" {
            continue;
        }
        let parent = gff_attr(cols[8], "Parent").unwrap_or("");
        let trans_id = parent.split(',').next().unwrap_or(parent);
        let Some(tx) = trans.get_mut(trans_id) else {
            continue;
        };
        let (start, end): (i64, i64) = (cols[3].parse()?, cols[4].parse()?);
        if feature == "exon" {
            tx.e_starts.push(start - 1);
            tx.e_ends.push(end);
        } else {
            tx.c_starts.push(start - 1);
            tx.c_ends.push(end);
        }
    }

    let mut out = tama_io::create_writer(output)?;
    for trans_id in &order {
        let tx = &trans[trans_id];
        let mut e_starts = tx.e_starts.clone();
        let mut e_ends = tx.e_ends.clone();
        e_starts.sort_unstable();
        e_ends.sort_unstable();
        let (mut blocks, mut rel_starts) = (String::new(), String::new());
        for k in 0..e_starts.len() {
            if k > 0 {
                blocks.push(',');
                rel_starts.push(',');
            }
            blocks.push_str(&(e_ends[k] - e_starts[k]).to_string());
            rel_starts.push_str(&(e_starts[k] - tx.t_start).to_string());
        }
        let (cds_start, cds_end) = if tx.c_starts.is_empty() {
            (0, 0)
        } else {
            let mut cs = tx.c_starts.clone();
            let mut ce = tx.c_ends.clone();
            cs.sort_unstable();
            ce.sort_unstable();
            (cs[0], *ce.last().unwrap())
        };
        writeln!(
            out,
            "{}\t{}\t{}\t{};{};{};{}\t40\t{}\t{cds_start}\t{cds_end}\t255,0,0\t{}\t{blocks}\t{rel_starts}",
            tx.chrom, tx.t_start, tx.t_end, tx.gene_id, trans_id, tx.trans_name, tx.trans_type,
            tx.strand, e_starts.len()
        )?;
    }
    Ok(())
}

/// Convert a Cupcake `collapsed.gff` to TAMA BED12. Ports
/// `tama_format_gff_to_bed12_cupcake`. Exons come from `exon` feature lines in
/// file order (which must be genomic-ascending).
fn gff2bed_cupcake(gff: &std::path::Path, output: &std::path::Path) -> anyhow::Result<()> {
    use std::io::BufRead;
    let mut gene_order: Vec<String> = Vec::new();
    let mut gene_trans: IndexMap<String, Vec<String>> = IndexMap::new();
    // (gene, trans) -> (chrom, strand, exons in order)
    type CupTx = (String, char, Vec<(i64, i64)>);
    let mut trans: IndexMap<(String, String), CupTx> = IndexMap::new();

    let reader = tama_io::open_reader(gff)?;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        let chrom = cols[0];
        let feature = cols[2];
        let e_start: i64 = cols[3].parse()?;
        let e_end: i64 = cols[4].parse()?;
        let strand = cols[6].chars().next().unwrap_or('+');

        let mut gene_id = None;
        let mut trans_id = None;
        for field in cols[8].split(';') {
            let field = field.trim();
            if field.is_empty() || !field.contains('"') {
                continue;
            }
            let id_code = field.split('"').nth(1).unwrap_or("");
            if field.contains("gene_id") {
                gene_id = Some(id_code.to_string());
            } else if field.contains("transcript_id") {
                trans_id = Some(id_code.to_string());
            }
        }
        let gene_id = gene_id.ok_or_else(|| anyhow::anyhow!("missing gene ID: {line}"))?;
        let trans_id = trans_id.ok_or_else(|| anyhow::anyhow!("missing transcript ID: {line}"))?;

        if !gene_trans.contains_key(&gene_id) {
            gene_order.push(gene_id.clone());
            gene_trans.insert(gene_id.clone(), Vec::new());
        }
        let key = (gene_id.clone(), trans_id.clone());
        if !trans.contains_key(&key) {
            gene_trans.get_mut(&gene_id).unwrap().push(trans_id.clone());
            trans.insert(key.clone(), (chrom.to_string(), strand, Vec::new()));
        }
        if feature == "exon" {
            trans.get_mut(&key).unwrap().2.push((e_start, e_end));
        }
    }

    let mut records = Vec::new();
    for gene in &gene_order {
        for trans_id in &gene_trans[gene] {
            let (chrom, strand, exons) = &trans[&(gene.clone(), trans_id.clone())];
            records.push((
                chrom.clone(),
                *strand,
                gene.clone(),
                trans_id.clone(),
                exons.clone(),
            ));
        }
    }
    write_gtf_bed(&records, output)
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
            let mut new = vec![
                ens_gene,
                ens_trans,
                tama_gene.to_string(),
                tama_trans.to_string(),
            ];
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
            let mut new = vec![
                ens_gene,
                ens_trans,
                tama_gene.to_string(),
                tama_trans.to_string(),
            ];
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
