# TAMA (Rust)

A Rust rewrite of [TAMA](https://github.com/GenomeRIK/tama) (Transcriptome
Annotation by Modular Algorithms), a toolkit for processing long-read
(Iso-Seq / Nanopore) transcriptome data.

The original project is a collection of Python 2 scripts. This rewrite ships a
single `tama` binary exposing every tool as a subcommand, with native Rust I/O
(no Python / BioPython / samtools runtime dependency).

## Building

```sh
cargo build --release
# binary at target/release/tama
```

## Layout

- `crates/tama-core` — domain model, TAMA BED12 dialect, CIGAR handling,
  sequence utilities, and (later) the collapse/merge algorithms.
- `crates/tama-io` — format I/O (gzip-aware readers, FASTA, SAM/BAM).
- `crates/tama` — the `tama` CLI (clap) dispatching to each subcommand.

## Subcommand map

| Original script | Subcommand |
|---|---|
| `tama_collapse.py` | `tama collapse` |
| `tama_merge.py` | `tama merge` |
| `tama_variant_caller.py` | `tama variants call` |
| `tama_degradation_signature.py` | `tama stats degradation` |
| `tama_find_model_changes.py` | `tama stats model-changes` |
| `tama_sampling_saturation_curve.py` | `tama stats saturation` |
| `tama_filter_primary_transcripts_orf.py` | `tama filter primary-orf` |
| `tama_remove_fragment_models.py` | `tama filter fragments` |
| `tama_remove_polya_models_levels.py` | `tama filter polya` |
| `tama_remove_single_read_models_levels.py` | `tama filter single-read` |
| `tama_convert_bed_gtf_ensembl_no_cds.py` | `tama format bed2gtf` |
| `tama_convert_bed_gtf_ensembl_orf_nmd.py` | `tama format bed2gtf-orf` |
| `tama_convert_nanopore_fastq_fasta.py` | `tama format fastq2fasta` |
| `tama_format_gff_to_bed12_*.py` | `tama format gff2bed --source <cupcake\|liftoff>` |
| `tama_format_gtf_to_bed12_*.py` | `tama format gtf2bed --source <ensembl\|ncbi\|stringtie>` |
| `tama_format_id_filter.py` | `tama format id-filter` |
| `tama_orf_seeker.py` | `tama orf seek` |
| `tama_bed_extract_cds.py` | `tama orf extract-cds` |
| `tama_cds_regions_bed_add.py` | `tama orf add-cds` |
| `tama_orf_blastp_parser.py` | `tama orf blastp-parse` |
| `tama_read_support_collapse_cluster.py` | `tama support collapse-cluster` |
| `tama_read_support_levels.py` | `tama support levels` |
| `tama_read_support_merge_collapse.py` | `tama support merge-collapse` |
| `tama_flnc_polya_cleanup.py` | `tama cleanup polya` |
| `tama_fasta_splitter.py` | `tama split fasta` |
| `tama_mapped_sam_splitter.py` | `tama split sam` |

## Status

Under active development. Implemented so far: CLI scaffold, TAMA BED12
read/write, CIGAR/coordinate handling. See the plan for the remaining phases.

## License

GPL-3.0-or-later, matching upstream TAMA.
