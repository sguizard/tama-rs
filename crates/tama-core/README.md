# tama-core

Core domain model and algorithms for [`tama-rs`](https://crates.io/crates/tama-rs),
a Rust rewrite of [TAMA](https://github.com/GenomeRIK/tama) (Transcriptome
Annotation by Modular Algorithms).

This crate holds the parts with no I/O: the TAMA BED12 dialect, CIGAR handling,
per-read error and variation calculation, poly-A detection, and the
collapse/merge grouping algorithms.

It is published so the algorithms can be reused; most users want the
[`tama-rs`](https://crates.io/crates/tama-rs) command-line tool instead.

## License

GPL-3.0-or-later, matching the original TAMA.
