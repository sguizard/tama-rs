//! `tama collapse` — collapse mapped long reads into transcript models.
//!
//! Ports `tama_collapse.py` (original run mode, capped path). Argument names
//! mirror the original short flags.

use std::io::Write;

use anyhow::{bail, Context};
use clap::Parser;
use indexmap::IndexMap;

use tama_core::collapse::{self, Cap, CollapseParams, CollapseResult, Ends, Merged, Transcript};
use tama_core::error_calc::{calc_error_rate, Variation};
use tama_core::gene::{gene_group, GeneMember};
use tama_core::metrics::{read_metrics, IdentMethod};
use tama_core::polya::{detect_polya, PolyA};
use tama_io::sam::{mapped_flag, read_sam};

const A_WINDOW: i64 = 20;
const A_PERC_THRESH: f64 = 70.0;

#[derive(Parser)]
pub struct Args {
    /// Sorted SAM file (or BAM with --bam). (`-s`)
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
    /// Collapse exon ends: `common_ends` or `longest_ends`. (`-e`)
    #[arg(short = 'e', long, default_value = "common_ends")]
    pub ends: String,
    /// Minimum coverage percent. (`-c`)
    #[arg(short = 'c', long, default_value_t = 99.0)]
    pub coverage: f64,
    /// Minimum identity percent. (`-i`)
    #[arg(short = 'i', long, default_value_t = 85.0)]
    pub identity: f64,
    /// Identity calculation method: `ident_cov` or `ident_map`. (`-icm`)
    #[arg(long = "icm", default_value = "ident_cov")]
    pub ident_method: String,
    /// 5' threshold. (`-a`)
    #[arg(short = 'a', long, default_value_t = 10)]
    pub five_prime: i64,
    /// Exon/splice-junction threshold. (`-m`)
    #[arg(short = 'm', long, default_value_t = 10)]
    pub exon_thresh: i64,
    /// 3' threshold. (`-z`)
    #[arg(short = 'z', long, default_value_t = 10)]
    pub three_prime: i64,
    /// Duplicate merge behaviour: `merge_dup` or `no_merge`. (`-d`)
    #[arg(short = 'd', long, default_value = "merge_dup")]
    pub dup: String,
    /// Splice-junction priority: `no_priority` or `sj_priority`. (`-sj`)
    #[arg(long = "sj", default_value = "no_priority")]
    pub sj_priority: String,
    /// Splice-junction error threshold (bp). (`-sjt`)
    #[arg(long = "sjt", default_value_t = 10)]
    pub sj_thresh: i64,
    /// Local density error threshold. (`-lde`)
    #[arg(long = "lde", default_value_t = 1000)]
    pub lde: i64,
    /// Treat input as BAM instead of SAM. (`-b`)
    #[arg(short = 'b', long)]
    pub bam: bool,
    /// Run mode: `original` or `low_mem`. (`-rm`)
    #[arg(long = "rm", default_value = "original")]
    pub run_mode: String,
    /// Variation coverage threshold (reads). (`-vc`)
    #[arg(long = "vc", default_value_t = 5)]
    pub var_coverage: i64,
}

/// Per-accepted-read state carried into collapsing/output.
struct ReadData {
    trans: Transcript,
    polya: PolyA,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    if args.bam {
        bail!("BAM input (-b) is not implemented yet; convert to SAM first");
    }
    let cap = match args.cap_flag.as_str() {
        "capped" => Cap::Capped,
        "no_cap" => Cap::NoCap,
        other => bail!("invalid -x cap flag: {other:?} (use capped or no_cap)"),
    };
    let ends = match args.ends.as_str() {
        "common_ends" => Ends::CommonEnds,
        "longest_ends" => Ends::LongestEnds,
        other => bail!("invalid -e ends flag: {other:?}"),
    };
    let ident_method = match args.ident_method.as_str() {
        "ident_cov" => IdentMethod::IdentCov,
        "ident_map" => IdentMethod::IdentMap,
        other => bail!("invalid -icm method: {other:?}"),
    };
    let params = CollapseParams {
        fiveprime_threshold: args.five_prime,
        threeprime_threshold: args.three_prime,
        exon_diff_threshold: args.exon_thresh,
        cap,
        ends,
        sj_priority: args.sj_priority == "sj_priority",
        merge_dup: args.dup != "no_merge",
    };

    let genome = tama_io::fasta::load_fasta(&args.fasta)
        .with_context(|| format!("loading genome {}", args.fasta.display()))?;
    let records =
        read_sam(&args.sam).with_context(|| format!("reading SAM {}", args.sam.display()))?;

    // ---- output writers ----
    let p = &args.prefix;
    let mut out_bed = tama_io::create_writer(format!("{p}.bed"))?;
    let mut out_read = tama_io::create_writer(format!("{p}_read.txt"))?;
    let mut out_trans_report = tama_io::create_writer(format!("{p}_trans_report.txt"))?;
    let mut out_trans_read = tama_io::create_writer(format!("{p}_trans_read.bed"))?;
    let mut out_polya = tama_io::create_writer(format!("{p}_polya.txt"))?;
    let mut out_strand = tama_io::create_writer(format!("{p}_strand_check.txt"))?;
    let mut out_report = tama_io::create_writer(format!("{p}_report.txt"))?;

    writeln!(out_read, "read_id\tmapped_flag\taccept_flag\tpercent_coverage\tpercent_identity\terror_line<h;s;i;d;m>\tlength\tcigar")?;
    writeln!(out_trans_report, "transcript_id\tnum_clusters\thigh_coverage\tlow_coverage\thigh_quality_percent\tlow_quality_percent\tstart_wobble_list\tend_wobble_list\tcollapse_sj_start_err\tcollapse_sj_end_err\tcollapse_error_nuc")?;
    writeln!(
        out_polya,
        "cluster_id\ttrans_id\tstrand\ta_percent\ta_count\tsequence"
    )?;
    writeln!(out_strand, "read_id\tscaff_name\tstart_pos\tcigar\tstrands")?;

    // ---- per-read pass: error/coverage/identity, acceptance, poly-A ----
    let mut variation = Variation::default();
    let mut accepted: IndexMap<String, ReadData> = IndexMap::new();
    let mut n_accepted = 0i64;
    let mut n_discarded = 0i64;

    for rec in &records {
        let mflag = mapped_flag(rec.flag);
        let strand = match mflag {
            "forward_strand" => '+',
            "reverse_strand" => '-',
            _ => {
                // unmapped / not_primary / chimeric
                writeln!(
                    out_read,
                    "{}\t{}\t{}\tNA\tNA\tNA\tNA\tNA",
                    rec.read_id, mflag, mflag
                )?;
                n_discarded += 1;
                continue;
            }
        };

        // strand-check vs XS tag
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
            n_discarded += 1;
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
        n_accepted += 1;

        let polya = detect_polya(genome_seq, strand, rec.start_pos, tc.end_pos, A_WINDOW);

        if accepted.contains_key(&rec.read_id) {
            bail!(
                "multi-mapped read {} — multimap handling not implemented yet",
                rec.read_id
            );
        }
        accepted.insert(
            rec.read_id.clone(),
            ReadData {
                trans: Transcript {
                    cluster_id: rec.read_id.clone(),
                    scaff_name: rec.scaffold.clone(),
                    strand,
                    start_pos: rec.start_pos,
                    end_pos: tc.end_pos,
                    exon_start_list: tc.exon_start_list,
                    exon_end_list: tc.exon_end_list,
                    sj_pre_error_list: er.sj_pre_error_list,
                    sj_post_error_list: er.sj_post_error_list,
                    percent_coverage: m.percent_coverage,
                    percent_identity: m.percent_identity,
                },
                polya,
            },
        );
    }

    // ---- position grouping (contiguous overlap on scaffold, SAM order) ----
    let groups = position_groups(&accepted)?;

    // ---- per group: gene grouping, collapsing, output ----
    let mut gene_count = 0i64;
    let mut report_trans_count = 0i64;

    for group in &groups {
        // split by strand
        let mut fwd: Vec<&str> = Vec::new();
        let mut rev: Vec<&str> = Vec::new();
        for id in group {
            match accepted[id.as_str()].trans.strand {
                '+' => fwd.push(id),
                '-' => rev.push(id),
                _ => {}
            }
        }
        let fwd_genes = make_genes(&fwd, &accepted);
        let rev_genes = make_genes(&rev, &accepted);

        // merge forward/reverse gene lists ordered by gene_start (forward first)
        let mut gene_entries: Vec<(i64, u8, Vec<String>)> = Vec::new();
        for g in fwd_genes {
            gene_entries.push((g.0, 0, g.1));
        }
        for g in rev_genes {
            gene_entries.push((g.0, 1, g.1));
        }
        gene_entries.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        for (_start, _rev, trans_ids) in gene_entries {
            gene_count += 1;
            let trans_objs: Vec<Transcript> = trans_ids
                .iter()
                .map(|id| accepted[id.as_str()].trans.clone())
                .collect();

            let match_groups = match cap {
                Cap::Capped => collapse::simplify_gene_capped(&trans_objs, &params),
                Cap::NoCap => tama_core::collapse_nocap::simplify_gene_nocap(&trans_objs, &params),
            };

            let mut merged_list: Vec<Merged> = Vec::new();
            for idx_group in &match_groups {
                let members: Vec<Transcript> =
                    idx_group.iter().map(|&k| trans_objs[k].clone()).collect();
                let cr = collapse::collapse_transcripts(&members, &params);
                merged_list.push(build_merged(&members, cr));
            }

            let sorted = collapse::sort_transcripts(merged_list, &params);
            let mut trans_count = 0i64;
            for mut merged in sorted {
                trans_count += 1;
                merged.trans_id = format!("G{gene_count}.{trans_count}");
                report_trans_count += 1;

                writeln!(out_bed, "{}", merged.format_bed_line())?;
                writeln!(
                    out_trans_report,
                    "{}",
                    format_trans_report(&merged, &accepted)
                )?;

                for cid in &merged.merged_trans {
                    let rd = &accepted[cid.as_str()];
                    writeln!(
                        out_trans_read,
                        "{}",
                        read_bed_line(&rd.trans, &merged.trans_id)
                    )?;

                    let a_percent = rd.polya.a_percent * 100.0;
                    if a_percent > A_PERC_THRESH {
                        writeln!(
                            out_polya,
                            "{}\t{}\t{}\t{}\t{}\t{}",
                            cid,
                            merged.trans_id,
                            rd.trans.strand,
                            round2(a_percent),
                            rd.polya.a_count,
                            rd.polya.downstream_seq
                        )?;
                    }
                }
            }
        }
    }

    // ---- report ----
    writeln!(out_report, "TAMA Collapse has run successfully!")?;
    writeln!(out_report, "Total Gene Count:\t{gene_count}")?;
    writeln!(out_report, "Total Transcript Count:\t{report_trans_count}")?;
    writeln!(out_report, "Total Accepted Reads:\t{n_accepted}")?;
    writeln!(out_report, "Total Discarded Reads:\t{n_discarded}")?;

    log::info!(
        "collapse done: {gene_count} genes, {report_trans_count} transcripts, {n_accepted} accepted, {n_discarded} discarded"
    );
    log::warn!("variant calling (_variants.txt/_varcov.txt) and _local_density_error.txt are not yet ported");
    Ok(())
}

/// Round to 2 decimal places, dropping trailing zeros like Python's `str(round(x,2))`.
fn round2(x: f64) -> String {
    let r = (x * 100.0).round() / 100.0;
    if r.fract() == 0.0 {
        format!("{:.1}", r)
    } else {
        let s = format!("{:.2}", r);
        s.trim_end_matches('0').to_string()
    }
}

/// Contiguous position grouping of accepted reads in SAM order.
fn position_groups(accepted: &IndexMap<String, ReadData>) -> anyhow::Result<Vec<Vec<String>>> {
    let mut groups: Vec<Vec<String>> = Vec::new();
    let mut this_scaffold = String::new();
    let mut group_start = 0i64;
    let mut group_end = 0i64;

    for (id, rd) in accepted {
        let t = &rd.trans;
        if this_scaffold.is_empty() {
            this_scaffold = t.scaff_name.clone();
            group_start = t.start_pos;
            group_end = t.end_pos;
            groups.push(vec![id.clone()]);
            continue;
        }
        if t.scaff_name == this_scaffold {
            if t.start_pos >= group_start && t.start_pos <= group_end {
                groups.last_mut().unwrap().push(id.clone());
                if t.end_pos > group_end {
                    group_end = t.end_pos;
                }
            } else if t.start_pos > group_end {
                group_start = t.start_pos;
                group_end = t.end_pos;
                groups.push(vec![id.clone()]);
            } else {
                bail!("SAM file not sorted at read {id}");
            }
        } else {
            this_scaffold = t.scaff_name.clone();
            group_start = t.start_pos;
            group_end = t.end_pos;
            groups.push(vec![id.clone()]);
        }
    }
    Ok(groups)
}

/// Gene-group a strand's reads, returning (gene_start, member ids) ordered by start.
fn make_genes(ids: &[&str], accepted: &IndexMap<String, ReadData>) -> Vec<(i64, Vec<String>)> {
    if ids.is_empty() {
        return Vec::new();
    }
    let members: Vec<GeneMember> = ids
        .iter()
        .map(|id| {
            let t = &accepted[*id].trans;
            GeneMember {
                id,
                exon_starts: &t.exon_start_list,
                exon_ends: &t.exon_end_list,
            }
        })
        .collect();
    gene_group(&members)
        .into_iter()
        .map(|g| (g.gene_start, g.trans_ids))
        .collect()
}

/// Build a `Merged` from members and a collapse result.
fn build_merged(members: &[Transcript], cr: CollapseResult) -> Merged {
    let num_exons = members.iter().map(|t| t.num_exons()).max().unwrap_or(0);
    Merged {
        trans_id: String::new(),
        scaff_name: members[0].scaff_name.clone(),
        strand: members[0].strand,
        start_pos: cr.collapse_start_list[0],
        end_pos: *cr.collapse_end_list.last().unwrap(),
        num_exons,
        collapse_start_list: cr.collapse_start_list,
        collapse_end_list: cr.collapse_end_list,
        start_wobble_list: cr.start_wobble_list,
        end_wobble_list: cr.end_wobble_list,
        collapse_sj_start_err_list: cr.collapse_sj_start_err_list,
        collapse_sj_end_err_list: cr.collapse_sj_end_err_list,
        collapse_start_error_nuc_list: cr.collapse_start_error_nuc_list,
        collapse_end_error_nuc_list: cr.collapse_end_error_nuc_list,
        merged_trans: members.iter().map(|t| t.cluster_id.clone()).collect(),
    }
}

/// Read-level BED line for a member transcript (mirrors `Transcript.format_bed_line`).
fn read_bed_line(t: &Transcript, final_trans_id: &str) -> String {
    let id_line = format!("{};{}", final_trans_id, t.cluster_id);
    let mut sizes = String::new();
    let mut starts = String::new();
    for k in 0..t.num_exons() {
        if k > 0 {
            sizes.push(',');
            starts.push(',');
        }
        sizes.push_str(&(t.exon_end_list[k] - t.exon_start_list[k]).to_string());
        starts.push_str(&(t.exon_start_list[k] - t.start_pos).to_string());
    }
    [
        t.scaff_name.clone(),
        (t.start_pos - 1).to_string(),
        (t.end_pos - 1).to_string(),
        id_line,
        "40".to_string(),
        t.strand.to_string(),
        (t.start_pos - 1).to_string(),
        (t.end_pos - 1).to_string(),
        "255,0,0".to_string(),
        t.num_exons().to_string(),
        sizes,
        starts,
    ]
    .join("\t")
}

/// Trans report line (mirrors `Merged.format_trans_report_line`).
fn format_trans_report(m: &Merged, accepted: &IndexMap<String, ReadData>) -> String {
    let mut hi_q = -1.0f64;
    let mut lo_q = -1.0f64;
    let mut hi_c = -1.0f64;
    let mut lo_c = -1.0f64;
    for cid in &m.merged_trans {
        let t = &accepted[cid.as_str()].trans;
        let q = t.percent_identity;
        let c = t.percent_coverage;
        if hi_q == -1.0 {
            hi_q = q;
            lo_q = q;
            hi_c = c;
            lo_c = c;
        } else {
            if hi_q < q {
                hi_q = q;
            }
            if lo_q > q {
                lo_q = q;
            }
            if hi_c < c {
                hi_c = c;
            }
            if lo_c > c {
                lo_c = c;
            }
        }
    }

    let join_i = |v: &[i64]| {
        v.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    };
    let (mut sj_start, mut sj_end) = (
        m.collapse_sj_start_err_list.clone(),
        m.collapse_sj_end_err_list.clone(),
    );
    if m.strand == '+' {
        sj_start.reverse();
        sj_end.reverse();
    }

    // collapse_error_nuc: drop the terminal "na" (5' start has no SJ) and, on the
    // plus strand, reverse into genomic 5'->3' order. Mirrors format_trans_report_line.
    let mut nuc = m.collapse_start_error_nuc_list.clone();
    if nuc.len() > 1 {
        if m.strand == '+' {
            if nuc.last().map(String::as_str) == Some("na") {
                nuc.pop();
                nuc.reverse();
            }
        } else if nuc.first().map(String::as_str) == Some("na") {
            nuc.remove(0);
        }
    }

    [
        m.trans_id.clone(),
        m.merged_trans.len().to_string(),
        round2(hi_c),
        round2(lo_c),
        round2(hi_q),
        round2(lo_q),
        join_i(&m.start_wobble_list),
        join_i(&m.end_wobble_list),
        join_i(&sj_start),
        join_i(&sj_end),
        nuc.join(";"),
    ]
    .join("\t")
}
