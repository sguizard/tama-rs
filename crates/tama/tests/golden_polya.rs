//! Golden test for `detect_polya` against the original TAMA `collapse_polya.txt`.
//! The reference run flags read `11_c69636/1/797` (minus strand) with a 20 bp
//! downstream window that is 75% A.

use std::path::PathBuf;

use tama_core::cigar::trans_coordinates;
use tama_core::polya::detect_polya;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..")
}

#[test]
fn polya_window_matches_golden() {
    let root = workspace_root();
    let genome_map =
        tama_io::fasta::load_fasta(root.join("test_data/test_genome.fa")).expect("genome");
    let sam = std::fs::read_to_string(root.join("test_data/gmap_test.sam")).expect("sam");

    let rec = sam
        .lines()
        .find(|l| l.starts_with("11_c69636/1/797\t"))
        .expect("read present");
    let f: Vec<&str> = rec.split('\t').collect();
    let scaff = f[2];
    let start_pos: i64 = f[3].parse().unwrap();
    let cigar = f[5];
    let genome = genome_map.get(scaff).unwrap().as_bytes();

    let tc = trans_coordinates(start_pos, cigar).unwrap();
    // SAM flag 16 => minus strand.
    let p = detect_polya(genome, '-', start_pos, tc.end_pos, 20);

    assert_eq!(p.downstream_seq, "AAGGAAAGAGGAAAAAAAAA");
    assert_eq!(p.a_count, 15);
    assert!((p.a_percent - 0.75).abs() < 1e-9);
}
