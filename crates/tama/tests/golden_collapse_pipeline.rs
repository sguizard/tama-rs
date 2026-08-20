//! Golden test for the full `tama collapse` pipeline (capped mode) against the
//! original TAMA outputs in `tests/golden/`.
//!
//! Asserts byte-identical output for `.bed`, `_read.txt`, `_polya.txt`,
//! `_strand_check.txt`, and `_trans_read.bed`. For `_trans_report.txt` the first
//! ten columns must match exactly; the final `collapse_error_nuc` column is
//! compared as a token set because the original's within-cell order is a Python-2
//! hash-order artifact.

use std::path::PathBuf;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..")
}

fn sorted_lines(s: &str) -> Vec<&str> {
    let mut v: Vec<&str> = s.lines().collect();
    v.sort_unstable();
    v
}

#[test]
fn collapse_pipeline_matches_golden() {
    let root = workspace_root();
    let out_dir = std::env::temp_dir().join(format!("tama_collapse_test_{}", std::process::id()));
    std::fs::create_dir_all(&out_dir).unwrap();
    let prefix = out_dir.join("collapse");

    let args = tama::cmd::collapse::Args {
        sam: root.join("test_data/gmap_test.sam"),
        fasta: root.join("test_data/test_genome.fa"),
        prefix: prefix.to_str().unwrap().to_string(),
        cap_flag: "capped".to_string(),
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

    let golden = root.join("tests/golden");
    let read = |p: PathBuf| std::fs::read_to_string(p).unwrap();

    // Exact, in-order match.
    assert_eq!(
        read(out_dir.join("collapse.bed")),
        read(golden.join("collapse.bed")),
        ".bed must match the original exactly (order + IDs)"
    );

    // Order-insensitive exact match.
    for name in ["collapse_read.txt", "collapse_polya.txt", "collapse_trans_read.bed"] {
        let mine = read(out_dir.join(name));
        let gold = read(golden.join(name));
        assert_eq!(sorted_lines(&mine), sorted_lines(&gold), "{name} must match");
    }

    // strand_check is header-only in this dataset.
    assert_eq!(
        read(out_dir.join("collapse_strand_check.txt")),
        read(golden.join("collapse_strand_check.txt"))
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
    assert_eq!(norm(&mine), norm(&gold), "trans_report columns/tokens must match");

    let _ = std::fs::remove_dir_all(&out_dir);
}
