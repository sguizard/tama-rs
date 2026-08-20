//! Golden tests for the Phase 4 standalone tools (format converters, splitters,
//! read support) against outputs from the original Python 2 scripts. Inputs live
//! in `test_data/`, golden outputs in `tests/golden_p4/`.
//!
//! These drive the CLI binary so the whole clap dispatch + tool runs end to end.

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

#[test]
fn bed2gtf_matches_golden() {
    let r = root();
    let out = std::env::temp_dir().join(format!("p4_bed2gtf_{}.gtf", std::process::id()));
    let status = tama()
        .args(["format", "bed2gtf"])
        .arg(r.join("tests/golden/collapse.bed"))
        .arg(&out)
        .status()
        .unwrap();
    assert!(status.success());
    assert_eq!(read(out.clone()), read(r.join("tests/golden_p4/bed2gtf.gtf")));
    let _ = std::fs::remove_file(out);
}

#[test]
fn gtf2bed_stringtie_matches_golden() {
    let r = root();
    let out = std::env::temp_dir().join(format!("p4_st_{}.bed", std::process::id()));
    let status = tama()
        .args(["format", "gtf2bed", "--source", "stringtie"])
        .arg(r.join("test_data/p4/stringtie.gtf"))
        .arg(&out)
        .status()
        .unwrap();
    assert!(status.success());
    assert_eq!(read(out.clone()), read(r.join("tests/golden_p4/gtf2bed_stringtie.bed")));
    let _ = std::fs::remove_file(out);
}

fn gtf2bed_case(source: &str, input: &str, golden: &str) {
    let r = root();
    let out = std::env::temp_dir().join(format!("p4_{source}_{}.bed", std::process::id()));
    let status = tama()
        .args(["format", "gtf2bed", "--source", source])
        .arg(r.join(input))
        .arg(&out)
        .status()
        .unwrap();
    assert!(status.success());
    assert_eq!(read(out.clone()), read(r.join(golden)), "gtf2bed {source}");
    let _ = std::fs::remove_file(out);
}

#[test]
fn gtf2bed_ensembl_matches_golden() {
    gtf2bed_case("ensembl", "test_data/p4/ensembl.gtf", "tests/golden_p4/gtf2bed_ensembl.bed");
}

#[test]
fn gtf2bed_ncbi_matches_golden() {
    gtf2bed_case("ncbi", "test_data/p4/ncbi.gtf", "tests/golden_p4/gtf2bed_ncbi.bed");
}

#[test]
fn gff2bed_cupcake_matches_golden() {
    let r = root();
    let out = std::env::temp_dir().join(format!("p4_cc_{}.bed", std::process::id()));
    let status = tama()
        .args(["format", "gff2bed", "--source", "cupcake"])
        .arg(r.join("test_data/p4/cupcake.gff"))
        .arg(&out)
        .status()
        .unwrap();
    assert!(status.success());
    assert_eq!(read(out.clone()), read(r.join("tests/golden_p4/gff2bed_cupcake.bed")));
    let _ = std::fs::remove_file(out);
}

#[test]
fn fastq2fasta_matches_golden() {
    let r = root();
    let out = std::env::temp_dir().join(format!("p4_fq_{}.fa", std::process::id()));
    let status = tama()
        .args(["format", "fastq2fasta"])
        .arg(r.join("test_data/p4/test.fastq"))
        .arg(&out)
        .status()
        .unwrap();
    assert!(status.success());
    assert_eq!(read(out.clone()), read(r.join("tests/golden_p4/fastq2fasta.fa")));
    let _ = std::fs::remove_file(out);
}

#[test]
fn id_filter_matches_golden() {
    let r = root();
    let out = std::env::temp_dir().join(format!("p4_idf_{}.bed", std::process::id()));
    let status = tama()
        .args(["format", "id-filter", "-b"])
        .arg(r.join("tests/golden/collapse.bed"))
        .arg("-o")
        .arg(&out)
        .status()
        .unwrap();
    assert!(status.success());
    assert_eq!(read(out.clone()), read(r.join("tests/golden_p4/id_filter.bed")));
    let _ = std::fs::remove_file(out);
}

#[test]
fn split_fasta_matches_golden() {
    let r = root();
    let prefix = std::env::temp_dir().join(format!("p4_splitfa_{}", std::process::id()));
    let status = tama()
        .args(["split", "fasta"])
        .arg(r.join("test_data/p4/multi.fa"))
        .arg(&prefix)
        .arg("2")
        .status()
        .unwrap();
    assert!(status.success());
    for i in 1..=2 {
        let mine = read(PathBuf::from(format!("{}_{i}.fa", prefix.display())));
        let gold = read(r.join(format!("tests/golden_p4/split_fasta_{i}.fa")));
        assert_eq!(mine, gold, "split fasta chunk {i}");
        let _ = std::fs::remove_file(format!("{}_{i}.fa", prefix.display()));
    }
}

#[test]
fn split_sam_matches_golden() {
    let r = root();
    let prefix = std::env::temp_dir().join(format!("p4_splitsam_{}", std::process::id()));
    let status = tama()
        .args(["split", "sam"])
        .arg(r.join("test_data/gmap_test.sam"))
        .arg("2")
        .arg(&prefix)
        .status()
        .unwrap();
    assert!(status.success());
    // single scaffold -> one output file
    let mine = read(PathBuf::from(format!("{}_1.sam", prefix.display())));
    let gold = read(r.join("tests/golden_p4/split_sam_1.sam"));
    assert_eq!(mine, gold);
    let _ = std::fs::remove_file(format!("{}_1.sam", prefix.display()));
}

#[test]
fn support_merge_collapse_matches_golden() {
    let r = root();
    let dir = std::env::temp_dir().join(format!("p4_mc_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let tr = r.join("tests/golden/collapse_trans_read.bed");

    // build a collapse support file (no-cluster mode) for the merge source.
    let support = dir.join("mc_support.txt");
    assert!(tama()
        .args(["support", "collapse-cluster"])
        .arg(&tr)
        .arg(&tr)
        .arg(&support)
        .status()
        .unwrap()
        .success());

    // filelist: filename<TAB>prefix<TAB>dir/
    let filelist = dir.join("filelist.txt");
    std::fs::write(&filelist, format!("mc_support.txt\ttestmerge\t{}/\n", dir.display())).unwrap();

    let out = dir.join("out.txt");
    assert!(tama()
        .args(["support", "merge-collapse"])
        .arg(r.join("tests/golden_merge/merged_merge.txt"))
        .arg(&filelist)
        .arg(&out)
        .status()
        .unwrap()
        .success());

    // cols 1-5 exact; source_trans_line / source_read_line compared as sets.
    let norm = |s: &str| -> Vec<String> {
        let mut rows: Vec<String> = s
            .lines()
            .skip(1)
            .map(|l| {
                let c: Vec<&str> = l.split('\t').collect();
                let head = c[..5.min(c.len())].join("\t");
                let mut st: Vec<&str> = c.get(5).map(|x| x.split(',').collect()).unwrap_or_default();
                st.sort_unstable();
                let mut rd: Vec<&str> =
                    c.get(6).map(|x| x.split([';', ',']).collect()).unwrap_or_default();
                rd.sort_unstable();
                format!("{head}\t{}\t{}", st.join(","), rd.join(","))
            })
            .collect();
        rows.sort_unstable();
        rows
    };
    assert_eq!(
        norm(&read(out)),
        norm(&read(r.join("tests/golden_p4/support_merge_collapse.txt")))
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn support_collapse_cluster_matches_golden() {
    // no-cluster mode: trans_read.bed as both collapse and cluster input.
    let r = root();
    let out = std::env::temp_dir().join(format!("p4_sup_{}.txt", std::process::id()));
    let tr = r.join("tests/golden/collapse_trans_read.bed");
    let status = tama()
        .args(["support", "collapse-cluster"])
        .arg(&tr)
        .arg(&tr)
        .arg(&out)
        .status()
        .unwrap();
    assert!(status.success());

    // Counts must match exactly; the cluster_line member order is a Python-2
    // dict-order artifact, so compare it as a set.
    let norm = |s: &str| -> Vec<String> {
        let mut rows: Vec<String> = s
            .lines()
            .map(|l| {
                let c: Vec<&str> = l.split('\t').collect();
                let head = c[..c.len().min(4)].join("\t");
                let mut clusters: Vec<&str> =
                    c.get(4).map(|x| x.split(';').collect()).unwrap_or_default();
                clusters.sort_unstable();
                format!("{head}\t{}", clusters.join(";"))
            })
            .collect();
        rows.sort_unstable();
        rows
    };
    assert_eq!(
        norm(&read(out.clone())),
        norm(&read(r.join("tests/golden_p4/support_collapse_cluster.txt")))
    );
    let _ = std::fs::remove_file(out);
}
