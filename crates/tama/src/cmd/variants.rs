//! `tama variants` — variant calling.

use std::io::Write;

use anyhow::{bail, Context};
use clap::{Args as ClapArgs, Subcommand};
use indexmap::{IndexMap, IndexSet};

use tama_core::error_calc::{calc_error_rate, Variation};
use tama_core::metrics::{read_metrics, IdentMethod};
use tama_io::sam::{mapped_flag, read_sam};

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Call variants from a sorted SAM against the genome. (tama_variant_caller)
    Call(CallArgs),
}

#[derive(clap::Parser)]
pub struct CallArgs {
    /// Sorted SAM file. (`-s`)
    #[arg(short = 's', long = "sam")]
    pub sam: std::path::PathBuf,
    /// Genome FASTA file. (`-f`)
    #[arg(short = 'f', long = "fasta")]
    pub fasta: std::path::PathBuf,
    /// Output prefix. (`-p`)
    #[arg(short = 'p', long = "prefix")]
    pub prefix: String,
    /// Capped flag: `capped` or `no_cap`. (`-x`)
    #[arg(short = 'x', long, default_value = "no_cap")]
    pub cap_flag: String,
    /// Minimum coverage percent. (`-c`)
    #[arg(short = 'c', long, default_value_t = 99.0)]
    pub coverage: f64,
    /// Minimum identity percent. (`-i`)
    #[arg(short = 'i', long, default_value_t = 85.0)]
    pub identity: f64,
    /// Identity method: `ident_cov` or `ident_map`. (`-icm`)
    #[arg(long = "icm", default_value = "ident_cov")]
    pub ident_method: String,
    /// Splice-junction error threshold (bp). (`-sjt`)
    #[arg(long = "sjt", default_value_t = 10)]
    pub sj_thresh: i64,
    /// Variation coverage threshold (reads). (`-vc`)
    #[arg(long = "vc", default_value_t = 5)]
    pub var_coverage: usize,
    /// Treat input as BAM instead of SAM. (`-b`)
    #[arg(short = 'b', long)]
    pub bam: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    match args.cmd {
        Cmd::Call(a) => call(a),
    }
}

/// One accepted read: its scaffold and exon coordinates for the coverage pass.
struct AcceptedRead {
    read_id: String,
    scaffold: String,
    exon_start_list: Vec<i64>,
    exon_end_list: Vec<i64>,
}

fn call(args: CallArgs) -> anyhow::Result<()> {
    if args.bam {
        bail!("BAM input (-b) is not implemented yet; convert to SAM first");
    }
    let ident_method = match args.ident_method.as_str() {
        "ident_cov" => IdentMethod::IdentCov,
        "ident_map" => IdentMethod::IdentMap,
        other => bail!("invalid -icm method: {other:?}"),
    };

    let genome = tama_io::fasta::load_fasta(&args.fasta)
        .with_context(|| format!("loading genome {}", args.fasta.display()))?;
    let records =
        read_sam(&args.sam).with_context(|| format!("reading SAM {}", args.sam.display()))?;

    let p = &args.prefix;
    let mut out_read = tama_io::create_writer(format!("{p}_read.txt"))?;
    let mut out_strand = tama_io::create_writer(format!("{p}_strand_check.txt"))?;
    let mut out_variant = tama_io::create_writer(format!("{p}_variants.txt"))?;
    let mut out_varcov = tama_io::create_writer(format!("{p}_varcov.txt"))?;

    writeln!(out_read, "read_id\tmapped_flag\taccept_flag\tpercent_coverage\tpercent_identity\terror_line<h;s;i;d;m>\tlength\tcigar")?;
    writeln!(out_strand, "read_id\tscaff_name\tstart_pos\tcigar\tstrands")?;
    writeln!(
        out_variant,
        "scaffold\tposition\ttype\tref_allele\talt_allele\tcount\tcov_count\tcluster_list"
    )?;
    writeln!(out_varcov, "positions\toverlap_clusters")?;

    let mut variation = Variation::default();
    let mut accepted: Vec<AcceptedRead> = Vec::new();
    let mut scaffold_order: Vec<String> = Vec::new();
    let mut seen_scaffold: IndexSet<String> = IndexSet::new();

    for rec in &records {
        let mflag = mapped_flag(rec.flag);
        let strand = match mflag {
            "forward_strand" => '+',
            "reverse_strand" => '-',
            _ => {
                writeln!(
                    out_read,
                    "{}\t{}\t{}\tNA\tNA\tNA\tNA\tNA",
                    rec.read_id, mflag, mflag
                )?;
                continue;
            }
        };
        if let Some(xs) = rec.xs_strand {
            if strand == '+' && xs == '-' {
                writeln!(
                    out_strand,
                    "{}\t{}\t{}\t{}\t+-",
                    rec.read_id, rec.scaffold, rec.start_pos, rec.cigar
                )?;
            } else if strand == '-' && xs == '+' {
                writeln!(
                    out_strand,
                    "{}\t{}\t{}\t{}\t-+",
                    rec.read_id, rec.scaffold, rec.start_pos, rec.cigar
                )?;
            }
        }

        let genome_seq = genome
            .get(&rec.scaffold)
            .with_context(|| format!("scaffold {} not in genome", rec.scaffold))?
            .as_bytes();

        // records variation + per-position coverage for every mapped read
        let er = calc_error_rate(
            rec.start_pos,
            &rec.cigar,
            rec.seq.as_bytes(),
            genome_seq,
            &rec.scaffold,
            &rec.read_id,
            args.sj_thresh,
            &mut variation,
        )?;
        let tc = tama_core::cigar::trans_coordinates(rec.start_pos, &rec.cigar)?;
        let m = read_metrics(rec.seq.len() as i64, &er, ident_method);

        if m.percent_coverage < args.coverage || m.percent_identity < args.identity {
            writeln!(
                out_read,
                "{}\t{}\tdiscarded\t{}\t{}\t{}\t{}\t{}",
                rec.read_id,
                mflag,
                round2(m.percent_coverage),
                round2(m.percent_identity),
                m.error_line,
                m.length,
                rec.cigar
            )?;
            continue;
        }
        writeln!(
            out_read,
            "{}\t{}\taccepted\t{}\t{}\t{}\t{}\t{}",
            rec.read_id,
            mflag,
            round2(m.percent_coverage),
            round2(m.percent_identity),
            m.error_line,
            m.length,
            rec.cigar
        )?;

        if seen_scaffold.insert(rec.scaffold.clone()) {
            scaffold_order.push(rec.scaffold.clone());
        }
        accepted.push(AcceptedRead {
            read_id: rec.read_id.clone(),
            scaffold: rec.scaffold.clone(),
            exon_start_list: tc.exon_start_list,
            exon_end_list: tc.exon_end_list,
        });
    }

    // coverage pass: accepted reads spanning a variant position count toward it
    for ar in &accepted {
        if let Some(cov) = variation.coverage.get_mut(&ar.scaffold) {
            for (i, &es) in ar.exon_start_list.iter().enumerate() {
                let ee = ar.exon_end_list[i];
                for this_coord in es..ee {
                    if let Some(set) = cov.get_mut(&this_coord) {
                        set.insert(ar.read_id.clone());
                    }
                }
            }
        }
    }

    // emit variants + varcov per scaffold (original mode)
    const VAR_TYPES: [char; 5] = ['H', 'S', 'M', 'I', 'D'];
    for scaffold in &scaffold_order {
        let Some(pos_map) = variation.dict.get(scaffold) else {
            continue;
        };
        let genome_seq = genome[scaffold].as_bytes();
        let cov_map = &variation.coverage[scaffold];

        // cov_group -> ordered positions
        let mut cov_group: IndexMap<String, IndexSet<String>> = IndexMap::new();

        for (&var_pos, type_map) in pos_map {
            let mut cov_reads: Vec<String> = cov_map
                .get(&var_pos)
                .map(|s| s.iter().cloned().collect())
                .unwrap_or_default();
            cov_reads.sort();
            let var_coverage = cov_reads.len();

            if var_pos > genome_seq.len() as i64 || var_pos < 0 {
                continue;
            }
            let default_ref = (genome_seq[var_pos as usize] as char).to_string();

            let mut accept_pos = false;
            for vt in VAR_TYPES {
                let Some(seq_map) = type_map.get(&vt) else {
                    continue;
                };
                let ref_allele = if vt == 'M' {
                    default_ref.clone()
                } else {
                    "NA".to_string()
                };
                for (alt_seq, reads) in seq_map {
                    let count = reads.len();
                    if count >= args.var_coverage {
                        accept_pos = true;
                        let read_line = reads.iter().cloned().collect::<Vec<_>>().join(",");
                        writeln!(
                            out_variant,
                            "{scaffold}\t{var_pos}\t{vt}\t{ref_allele}\t{alt_seq}\t{count}\t{var_coverage}\t{read_line}"
                        )?;
                    }
                }
            }
            if accept_pos {
                let cov_line = cov_reads.join(",");
                cov_group
                    .entry(cov_line)
                    .or_default()
                    .insert(format!("{scaffold}_{var_pos}"));
            }
        }

        for (cov_line, positions) in &cov_group {
            let mut pos: Vec<String> = positions.iter().cloned().collect();
            pos.sort();
            writeln!(out_varcov, "{}\t{cov_line}", pos.join(","))?;
        }
    }

    log::info!("variant caller done: {} accepted reads", accepted.len());
    Ok(())
}

fn round2(x: f64) -> String {
    let r = (x * 100.0).round() / 100.0;
    if r.fract() == 0.0 {
        format!("{r:.1}")
    } else {
        format!("{r:.2}").trim_end_matches('0').to_string()
    }
}
