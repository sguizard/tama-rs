//! Golden test for `tama variants call` against the original Python 2
//! tama_variant_caller on `test_data/gmap_test.sam` (capped).
//!
//! `_read.txt`, `_strand_check.txt`, and `_varcov.txt` must match byte-for-byte;
//! `_variants.txt` matches on all columns with the `cluster_list` read order
//! compared as a set (that column reflects Python-2 dict iteration order).

use std::path::PathBuf;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..")
}

#[test]
fn variants_call_matches_golden() {
    let r = root();
    let dir = std::env::temp_dir().join(format!("vc_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let prefix = dir.join("vc");

    let status = std::process::Command::new(env!("CARGO_BIN_EXE_tama"))
        .args(["variants", "call", "-s"])
        .arg(r.join("test_data/gmap_test.sam"))
        .arg("-f")
        .arg(r.join("test_data/test_genome.fa"))
        .arg("-p")
        .arg(&prefix)
        .args(["-x", "capped"])
        .status()
        .unwrap();
    assert!(status.success());

    let read = |p: String| std::fs::read_to_string(p).unwrap();
    let p = prefix.display();

    for (suf, gold) in [
        ("_read.txt", "read.txt"),
        ("_strand_check.txt", "strand_check.txt"),
        ("_varcov.txt", "varcov.txt"),
    ] {
        assert_eq!(
            read(format!("{p}{suf}")),
            read(r.join(format!("tests/golden_variants/{gold}")).to_str().unwrap().to_string()),
            "{suf}"
        );
    }

    // variants: columns 1-7 exact, cluster_list (col 8) as a read set.
    let norm = |s: &str| -> Vec<String> {
        let mut rows: Vec<String> = s
            .lines()
            .filter(|l| !l.starts_with("scaffold"))
            .map(|l| {
                let c: Vec<&str> = l.split('\t').collect();
                let head = c[..7.min(c.len())].join("\t");
                let mut reads: Vec<&str> = c.get(7).map(|x| x.split(',').collect()).unwrap_or_default();
                reads.sort_unstable();
                format!("{head}\t{}", reads.join(","))
            })
            .collect();
        rows.sort_unstable();
        rows
    };
    assert_eq!(
        norm(&read(format!("{p}_variants.txt"))),
        norm(&read(r.join("tests/golden_variants/variants.txt").to_str().unwrap().to_string()))
    );
    let _ = std::fs::remove_dir_all(&dir);
}
