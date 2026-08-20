//! Golden test: validate the ported per-read error/coverage/identity pipeline
//! against reference output produced by the original Python 2 `tama_collapse.py`
//! run on `test_data/gmap_test.sam` (see `tests/golden/collapse_read.txt`).
//!
//! For every mapped read we recompute h/s/i/d/mismatch counts, coverage,
//! identity, length and the error line using `tama-core`, and assert they match
//! the golden values. This exercises `cigar`, `error_calc`, and `metrics`.

use std::path::PathBuf;

use tama_core::error_calc::{calc_error_rate, Variation};
use tama_core::metrics::{read_metrics, IdentMethod};

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
}

/// Map a SAM flag to TAMA's mapped-flag category.
fn mapped_flag(flag: u32) -> &'static str {
    match flag {
        0 => "forward_strand",
        16 => "reverse_strand",
        4 => "unmapped",
        256 | 272 => "not_primary",
        2048 | 2064 => "chimeric",
        _ => "unknown",
    }
}

fn round2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

#[test]
fn per_read_metrics_match_golden() {
    let root = workspace_root();
    let genome_map =
        tama_io::fasta::load_fasta(root.join("test_data/test_genome.fa")).expect("load genome");

    let sam = std::fs::read_to_string(root.join("test_data/gmap_test.sam")).expect("read sam");
    let golden =
        std::fs::read_to_string(root.join("tests/golden/collapse_read.txt")).expect("read golden");

    // Golden data lines in SAM order (skip header).
    let golden_lines: Vec<&str> = golden.lines().skip(1).collect();

    let sam_records: Vec<&str> = sam
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('@'))
        .collect();

    assert_eq!(
        sam_records.len(),
        golden_lines.len(),
        "record count mismatch between SAM and golden read.txt"
    );

    let mut checked = 0;
    for (rec, gline) in sam_records.iter().zip(&golden_lines) {
        let f: Vec<&str> = rec.split('\t').collect();
        let read_id = f[0];
        let flag: u32 = f[1].parse().unwrap();
        let scaff = f[2];
        let start_pos: i64 = f[3].parse().unwrap();
        let cigar = f[5];
        let read_seq = f[9];

        let g: Vec<&str> = gline.split('\t').collect();
        assert_eq!(g[0], read_id, "golden/SAM ordering diverged");

        let mflag = mapped_flag(flag);
        if mflag != "forward_strand" && mflag != "reverse_strand" {
            // Unmapped/not-primary/chimeric: golden records NA; just check the flag.
            assert_eq!(g[1], mflag, "mapped_flag mismatch for {read_id}");
            continue;
        }

        let genome = genome_map
            .get(scaff)
            .unwrap_or_else(|| panic!("scaffold {scaff} missing from genome"))
            .as_bytes();

        let mut var = Variation::default();
        let er = calc_error_rate(
            start_pos,
            cigar,
            read_seq.as_bytes(),
            genome,
            scaff,
            read_id,
            10,
            &mut var,
        )
        .expect("calc_error_rate");

        let m = read_metrics(read_seq.len() as i64, &er, IdentMethod::IdentCov);

        // error_line (integer counts) must match exactly.
        assert_eq!(m.error_line, g[5], "error_line mismatch for {read_id}");
        // length must match exactly.
        assert_eq!(m.length.to_string(), g[6], "length mismatch for {read_id}");

        // coverage / identity: compare rounded to 2 dp within tolerance.
        let gcov: f64 = g[3].parse().unwrap();
        let gident: f64 = g[4].parse().unwrap();
        assert!(
            (round2(m.percent_coverage) - gcov).abs() < 0.011,
            "coverage mismatch for {read_id}: rust={} golden={}",
            m.percent_coverage,
            gcov
        );
        assert!(
            (round2(m.percent_identity) - gident).abs() < 0.011,
            "identity mismatch for {read_id}: rust={} golden={}",
            m.percent_identity,
            gident
        );

        checked += 1;
    }

    assert!(checked > 50, "expected to check many reads, only did {checked}");
}
