//! Golden tests for Phase 5 tools (filters, stats, cleanup) and the read-support
//! `levels` tool, against the original Python 2 scripts.

use std::path::PathBuf;
use std::process::Command;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..")
}
fn tama() -> Command {
    Command::new(env!("CARGO_BIN_EXE_tama"))
}
fn read(p: PathBuf) -> String {
    std::fs::read_to_string(p).unwrap()
}

/// Normalize a read-support levels file: exact columns 1-5, `support_line`
/// (col 6) compared as per-source read sets (member order is a py2 artifact).
fn norm_levels(s: &str) -> Vec<String> {
    let mut rows: Vec<String> = s
        .lines()
        .filter(|l| !l.starts_with("merge_gene_id"))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            let head = c[..5.min(c.len())].join("\t");
            let mut groups: Vec<String> = c
                .get(5)
                .map(|x| {
                    x.split(';')
                        .map(|g| {
                            let (src, reads) = g.split_once(':').unwrap_or((g, ""));
                            let mut r: Vec<&str> = reads.split(',').collect();
                            r.sort_unstable();
                            format!("{src}:{}", r.join(","))
                        })
                        .collect()
                })
                .unwrap_or_default();
            groups.sort_unstable();
            format!("{head}\t{}", groups.join(";"))
        })
        .collect();
    rows.sort_unstable();
    rows
}

fn run_levels(merge: &str) -> String {
    let r = root();
    let dir = std::env::temp_dir().join(format!("p5_lv_{}_{}", merge.len(), std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let tr = r.join("tests/golden/collapse_trans_read.bed");
    let filelist = dir.join("filelist.txt");
    std::fs::write(&filelist, format!("testmerge\t{}\ttrans_read\n", tr.display())).unwrap();
    let prefix = dir.join("out");
    let merge_arg = if merge == "no_merge" {
        "no_merge".to_string()
    } else {
        r.join(merge).to_str().unwrap().to_string()
    };
    assert!(tama()
        .args(["support", "levels", "-f"])
        .arg(&filelist)
        .args(["-m", &merge_arg, "-o"])
        .arg(&prefix)
        .status()
        .unwrap()
        .success());
    let out = read(PathBuf::from(format!("{}_read_support.txt", prefix.display())));
    let _ = std::fs::remove_dir_all(&dir);
    out
}

#[test]
fn support_levels_no_merge_matches_golden() {
    let mine = run_levels("no_merge");
    let gold = read(root().join("tests/golden_p5/levels_no_merge.txt"));
    assert_eq!(norm_levels(&mine), norm_levels(&gold));
}

#[test]
fn support_levels_merge_matches_golden() {
    let mine = run_levels("tests/golden_merge/merged_merge.txt");
    let gold = read(root().join("tests/golden_p5/levels_merge.txt"));
    assert_eq!(norm_levels(&mine), norm_levels(&gold));
}

#[test]
fn filter_single_read_matches_golden() {
    let r = root();
    let dir = std::env::temp_dir().join(format!("p5_sr_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let prefix = dir.join("out");
    assert!(tama()
        .args(["filter", "single-read", "-b"])
        .arg(r.join("tests/golden/collapse.bed"))
        .arg("-r")
        .arg(r.join("tests/golden_p5/levels_no_merge.txt"))
        .arg("-o")
        .arg(&prefix)
        .status()
        .unwrap()
        .success());
    let p = prefix.display();
    assert_eq!(read(PathBuf::from(format!("{p}.bed"))), read(r.join("tests/golden_p5/single_read.bed")));
    assert_eq!(
        read(PathBuf::from(format!("{p}_singleton_report.txt"))),
        read(r.join("tests/golden_p5/single_read_report.txt"))
    );
    assert_eq!(
        read(PathBuf::from(format!("{p}_singleton.bed"))),
        read(r.join("tests/golden_p5/single_read_singleton.bed"))
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn filter_fragments_matches_golden() {
    let r = root();
    let dir = std::env::temp_dir().join(format!("p5_frag_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let prefix = dir.join("out");
    assert!(tama()
        .args(["filter", "fragments", "-f"])
        .arg(r.join("tests/golden_nocap/collapse.bed"))
        .arg("-o")
        .arg(&prefix)
        .status()
        .unwrap()
        .success());
    let p = prefix.display();
    assert_eq!(read(PathBuf::from(format!("{p}.bed"))), read(r.join("tests/golden_p5/fragments.bed")));
    assert_eq!(
        read(PathBuf::from(format!("{p}_discarded.txt"))),
        read(r.join("tests/golden_p5/fragments_discarded.txt"))
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn stats_degradation_matches_golden() {
    let r = root();
    let out = std::env::temp_dir().join(format!("p5_deg_{}.txt", std::process::id()));
    assert!(tama()
        .args(["stats", "degradation", "-c"])
        .arg(r.join("tests/golden/collapse_trans_read.bed"))
        .arg("--nc")
        .arg(r.join("tests/golden_nocap/collapse_trans_read.bed"))
        .arg("-o")
        .arg(&out)
        .status()
        .unwrap()
        .success());
    assert_eq!(read(out.clone()), read(r.join("tests/golden_p5/degradation.txt")));
    let _ = std::fs::remove_file(out);
}

#[test]
fn stats_saturation_structure() {
    // The saturation curve is read-order dependent (py2 dict order), so only the
    // structure (header, row count, read_count column) is deterministic.
    let r = root();
    let dir = std::env::temp_dir().join(format!("p5_sat_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let out = dir.join("sat.txt");
    assert!(tama()
        .args(["stats", "saturation", "-r"])
        .arg(r.join("tests/golden_p5/levels_no_merge.txt"))
        .args(["-b", "10", "-o"])
        .arg(&out)
        .status()
        .unwrap()
        .success());
    let content = read(out);
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines[0], "read_count\tgene_count");
    // read_count column is 10,20,30,... in order
    for (i, l) in lines[1..].iter().enumerate() {
        let rc: usize = l.split('\t').next().unwrap().parse().unwrap();
        assert_eq!(rc, (i + 1) * 10);
    }
    let _ = std::fs::remove_dir_all(&dir);
}
