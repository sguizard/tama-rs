//! Golden tests for `tama merge` with **no_cap** and **mixed** (capped + no_cap)
//! sources, against the original TAMA outputs in `tests/golden_merge_nocap/`.
//!
//! `.bed`, `_merge.txt`, and `_gene_report.txt` must match byte-for-byte
//! (order-insensitive) — these are the merged annotation and its provenance.
//!
//! `_trans_report.txt`: the deterministic columns (transcript_id, num_clusters,
//! sources, start/end wobble) must match exactly, and `all_source_trans` is
//! compared as a member set. The per-exon support columns are not asserted on:
//! when a coordinate is contended by a capped vs a no_cap member, which member
//! "wins" the support cell follows Python-2 dict iteration order and is not
//! reproducible (the same class of artifact documented for capped merge).

use std::path::{Path, PathBuf};

use tama::cmd::opts::{Dup, EndsOpt};

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
}

fn sorted_lines(s: &str) -> Vec<&str> {
    let mut v: Vec<&str> = s.lines().collect();
    v.sort_unstable();
    v
}

fn run_merge(filelist: &str, prefix: &Path) {
    let root = workspace_root();
    let args = tama::cmd::merge::Args {
        filelist: root.join(filelist),
        prefix: prefix.to_str().unwrap().to_string(),
        ends: EndsOpt::CommonEnds,
        five_prime: 20,
        exon_thresh: 10,
        three_prime: 20,
        dup: Dup::MergeDup,
        source_id: None,
        cds_source: None,
    };
    tama::cmd::merge::run(args).expect("merge run");
}

/// Keep the deterministic trans_report columns (transcript_id, num_clusters,
/// sources, start/end wobble) exactly, and reduce `all_source_trans` (col 7) to a
/// sorted member set. The per-exon support columns (5, 6) are dropped — their
/// contended-coordinate winner is a Python-2 dict-order artifact.
fn norm_trans_report(s: &str) -> Vec<String> {
    let mut out: Vec<String> = s
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            let mut members: Vec<&str> =
                f.get(7).map(|c| c.split(',').collect()).unwrap_or_default();
            members.sort_unstable();
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                f[0], // transcript_id
                f[1], // num_clusters
                f[2], // sources
                f[3], // start_wobble_list
                f[4], // end_wobble_list
                members.join(",")
            )
        })
        .collect();
    out.sort();
    out
}

fn check_case(case: &str, filelist: &str) {
    let root = workspace_root();
    let out_dir = std::env::temp_dir().join(format!("tama_merge_{case}_{}", std::process::id()));
    std::fs::create_dir_all(&out_dir).unwrap();
    let prefix = out_dir.join("m");
    run_merge(filelist, &prefix);

    let golden = root.join("tests/golden_merge_nocap");
    let read = |p: PathBuf| std::fs::read_to_string(p).unwrap();

    for (out, gold) in [
        ("m.bed", format!("{case}.bed")),
        ("m_merge.txt", format!("{case}_merge.txt")),
        ("m_gene_report.txt", format!("{case}_gene_report.txt")),
    ] {
        let mine = read(out_dir.join(out));
        let want = read(golden.join(&gold));
        assert_eq!(
            sorted_lines(&mine),
            sorted_lines(&want),
            "{gold} must match"
        );
    }

    let mine = read(out_dir.join("m_trans_report.txt"));
    let want = read(golden.join(format!("{case}_trans_report.txt")));
    assert_eq!(
        norm_trans_report(&mine),
        norm_trans_report(&want),
        "{case}_trans_report.txt must match (support/member cells as sets)"
    );
}

#[test]
fn merge_all_nocap_matches_golden() {
    check_case("nocap", "test_data/merge_nocap/fl_nocap.txt");
}

#[test]
fn merge_mixed_matches_golden() {
    check_case("mixed", "test_data/merge_nocap/fl_mixed.txt");
}
