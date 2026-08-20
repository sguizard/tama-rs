//! Golden test for the full `tama collapse` pipeline against the original TAMA
//! outputs, for both the capped (`tests/golden/`) and no-cap
//! (`tests/golden_nocap/`) modes.
//!
//! Asserts byte-identical output for `.bed` (in order), `_read.txt`,
//! `_polya.txt`, `_strand_check.txt`, and `_trans_read.bed`. For
//! `_trans_report.txt` the first ten columns must match exactly; the final
//! `collapse_error_nuc` column is compared as a token set because the original's
//! within-cell order is a Python-2 hash-order artifact.

use std::path::PathBuf;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..")
}

fn sorted_lines(s: &str) -> Vec<&str> {
    let mut v: Vec<&str> = s.lines().collect();
    v.sort_unstable();
    v
}

fn run_and_compare(cap_flag: &str, golden_subdir: &str) {
    let root = workspace_root();
    let out_dir = std::env::temp_dir().join(format!(
        "tama_collapse_test_{}_{}",
        cap_flag,
        std::process::id()
    ));
    std::fs::create_dir_all(&out_dir).unwrap();
    let prefix = out_dir.join("collapse");

    let args = tama::cmd::collapse::Args {
        sam: root.join("test_data/gmap_test.sam"),
        fasta: root.join("test_data/test_genome.fa"),
        prefix: prefix.to_str().unwrap().to_string(),
        cap_flag: cap_flag.to_string(),
        ends: "common_ends".to_string(),
        coverage: 99.0,
        identity: 85.0,
        ident_method: "ident_cov".to_string(),
        five_prime: 10,
        exon_thresh: 10,
        three_prime: 10,
        dup: "merge_dup".to_string(),
        sj_priority: "no_priority".to_string(),
        sj_thresh: 10,
        lde: 1000,
        bam: false,
        run_mode: "original".to_string(),
        var_coverage: 5,
    };
    tama::cmd::collapse::run(args).expect("collapse run");

    let golden = root.join(golden_subdir);
    let read = |p: PathBuf| std::fs::read_to_string(p).unwrap();

    assert_eq!(
        read(out_dir.join("collapse.bed")),
        read(golden.join("collapse.bed")),
        "[{cap_flag}] .bed must match the original exactly (order + IDs)"
    );

    for name in ["collapse_read.txt", "collapse_polya.txt", "collapse_trans_read.bed"] {
        let mine = read(out_dir.join(name));
        let gold = read(golden.join(name));
        assert_eq!(sorted_lines(&mine), sorted_lines(&gold), "[{cap_flag}] {name}");
    }

    assert_eq!(
        read(out_dir.join("collapse_strand_check.txt")),
        read(golden.join("collapse_strand_check.txt")),
        "[{cap_flag}] strand_check"
    );

    // trans_report: first 10 columns exact; last column as a token set.
    let mine = read(out_dir.join("collapse_trans_report.txt"));
    let gold = read(golden.join("collapse_trans_report.txt"));
    let norm = |s: &str| -> Vec<String> {
        let mut rows: Vec<String> = s
            .lines()
            .map(|l| {
                let cols: Vec<&str> = l.split('\t').collect();
                let head = cols[..cols.len().min(10)].join("\t");
                let mut toks: Vec<&str> = cols
                    .get(10)
                    .map(|c| c.split([';', '-']).collect())
                    .unwrap_or_default();
                toks.sort_unstable();
                format!("{head}\t{}", toks.join("|"))
            })
            .collect();
        rows.sort_unstable();
        rows
    };
    assert_eq!(norm(&mine), norm(&gold), "[{cap_flag}] trans_report");

    let _ = std::fs::remove_dir_all(&out_dir);
}

#[test]
fn collapse_capped_matches_golden() {
    run_and_compare("capped", "tests/golden");
}

#[test]
fn collapse_nocap_matches_golden() {
    run_and_compare("no_cap", "tests/golden_nocap");
}
