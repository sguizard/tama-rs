//! Golden test for `tama merge` (capped, single source) against the original
//! TAMA outputs in `tests/golden_merge/`.
//!
//! `.bed`, `_merge.txt`, and `_gene_report.txt` must match byte-for-byte
//! (order-insensitive). `_trans_report.txt` matches on all columns except the
//! member order within the `all_source_trans` column, which reflects Python-2
//! dict insertion order; it is compared as a member set.

use std::path::PathBuf;

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

#[test]
fn merge_capped_matches_golden() {
    let root = workspace_root();
    let out_dir = std::env::temp_dir().join(format!("tama_merge_test_{}", std::process::id()));
    std::fs::create_dir_all(&out_dir).unwrap();
    let prefix = out_dir.join("merged");

    let args = tama::cmd::merge::Args {
        filelist: root.join("test_data/merge/filelist_merge.txt"),
        prefix: prefix.to_str().unwrap().to_string(),
        ends: "common_ends".to_string(),
        five_prime: 20,
        exon_thresh: 10,
        three_prime: 20,
        dup: "no_merge".to_string(),
        source_id: None,
        cds_source: None,
    };
    tama::cmd::merge::run(args).expect("merge run");

    let golden = root.join("tests/golden_merge");
    let read = |p: PathBuf| std::fs::read_to_string(p).unwrap();

    for name in ["merged.bed", "merged_merge.txt", "merged_gene_report.txt"] {
        let mine = read(out_dir.join(name));
        let gold = read(golden.join(name));
        assert_eq!(
            sorted_lines(&mine),
            sorted_lines(&gold),
            "{name} must match"
        );
    }

    // trans_report: all columns exact except `all_source_trans` (col 8) which is
    // compared as a member set.
    let mine = read(out_dir.join("merged_trans_report.txt"));
    let gold = read(golden.join("merged_trans_report.txt"));
    let norm = |s: &str| -> Vec<String> {
        let mut rows: Vec<String> = s
            .lines()
            .map(|l| {
                let cols: Vec<&str> = l.split('\t').collect();
                let head = cols[..cols.len().min(7)].join("\t");
                let mut members: Vec<&str> = cols
                    .get(7)
                    .map(|c| c.split(',').collect())
                    .unwrap_or_default();
                members.sort_unstable();
                format!("{head}\t{}", members.join(","))
            })
            .collect();
        rows.sort_unstable();
        rows
    };
    assert_eq!(norm(&mine), norm(&gold), "trans_report must match");

    let _ = std::fs::remove_dir_all(&out_dir);
}
