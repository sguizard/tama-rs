# Golden reference outputs

These files are the output of the **original Python 2 `tama_collapse.py`** run on
`test_data/`, used as ground truth for the Rust integration tests in
`crates/tama/tests/`.

They are checked in so the tests run without a Python 2 environment. Regenerate
them with `./regenerate.sh` (requires conda/mamba to build a Python 2.7 +
BioPython env).

## Files

Produced by: `tama_collapse -s gmap_test.sam -f test_genome.fa -p collapse -x capped`

| File | Validated by |
|---|---|
| `collapse_read.txt` | `golden_collapse_read.rs` — per-read h/s/i/d/mismatch, coverage, identity, length |
| `collapse_polya.txt` | `golden_polya.rs` — poly-A window detection |
| `collapse.bed`, `collapse_trans_report.txt`, `collapse_trans_read.bed` | (future) collapse pipeline output |
| `collapse_variants.txt`, `collapse_varcov.txt` | (future) variant calling |
| `collapse_local_density_error.txt`, `collapse_report.txt`, `collapse_strand_check.txt` | (future) |
