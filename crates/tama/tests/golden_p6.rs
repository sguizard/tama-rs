//! Golden tests for Phase 6 tools (ORF/NMD) against the original Python 2 scripts.

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
fn orf_seek_matches_golden() {
    let r = root();
    let out = std::env::temp_dir().join(format!("p6_orf_{}.fa", std::process::id()));
    assert!(tama()
        .args(["orf", "seek", "-f"])
        .arg(r.join("test_data/p6/trans.fa"))
        .arg("-o")
        .arg(&out)
        .status()
        .unwrap()
        .success());
    assert_eq!(read(out.clone()), read(r.join("tests/golden_p6/orf_seek.fa")));
    let _ = std::fs::remove_file(out);
}
