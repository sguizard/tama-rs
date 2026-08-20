# TAMA (Rust)

A Rust rewrite of [TAMA](https://github.com/GenomeRIK/tama) (Transcriptome
Annotation by Modular Algorithms), a toolkit for processing long-read
(Iso-Seq / Nanopore) transcriptome data.

The original project is a collection of Python 2 scripts. This rewrite ships a
single `tama` binary exposing **every** tool as a subcommand, with native Rust
I/O and **no Python / BioPython / samtools runtime dependency**.

## Building

```sh
cargo build --release
# binary at target/release/tama
```

Requires a stable Rust toolchain (built and tested with 1.98).

## Usage

`tama <group> <tool>`. Every tool mirrors the original's short flags. Examples:

```sh
# collapse mapped long reads into transcript models
tama collapse -s sorted.sam -f genome.fa -p out -x capped

# merge annotations across sources
tama merge -f filelist.txt -p merged

# call variants from a sorted SAM
tama variants call -s sorted.sam -f genome.fa -p variants -x capped

# convert a TAMA BED to GTF
tama format bed2gtf annotation.bed annotation.gtf

# find ORFs, add CDS, and emit an ORF-annotated GTF
tama orf seek -f transcripts.fa -o orfs.faa
tama orf add-cds -p blastp_parse.txt -a annotation.bed -f transcripts.fa -o orf_nmd.bed
tama format bed2gtf-orf orf_nmd.bed orf_nmd.gtf
```

Run `tama --help`, `tama <group> --help`, or `tama <group> <tool> --help` for
the full flag list of any tool.

## Layout

- `crates/tama-core` — domain model, TAMA BED12 dialect, CIGAR handling, per-read
  error/variation calculation, and the collapse/merge grouping algorithms.
- `crates/tama-io` — format I/O (gzip-aware readers/writers, FASTA, SAM).
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

`gff2bed`/`gtf2bed` fold the per-source scripts into one subcommand with a
`--source` flag.

## Status

All 34 original TAMA tools are implemented.

Each tool is validated against the original Python 2 script by running both on
the same input and diffing the output ("golden" tests). Outputs are
**byte-identical** to the original, with one class of exception: a few
*diagnostic report columns* list the same tokens in a different order, because
the original derives that order from Python-2 dict iteration (hash order), which
is not reproducible and is not semantically meaningful. These are compared as
sets in the tests and are noted in code. Examples: `collapse` trans-report
`collapse_error_nuc`, `merge` `all_source_trans`, the read lists in `support`
outputs, the `variants` `cluster_list`, and `stats model-changes` diff cells.

Not yet ported (secondary collapse features): BAM input (`-b`), multimap
handling, the `_local_density_error.txt` output, and the `-rm low_mem` run mode.
Collapse's own `_variants.txt`/`_varcov.txt` are not written by `collapse`
itself yet, but the identical logic is available via `tama variants call`.

## Testing

```sh
cargo test          # unit + golden integration tests
cargo clippy --all-targets
cargo fmt --check
```

The golden reference outputs live under `tests/golden*`. To regenerate them from
the original Python 2 tools you need a Python 2.7 + BioPython environment; see
`tests/golden/regenerate.sh` (it builds one with conda/mamba).

## License

GPL-3.0-or-later, matching upstream TAMA.
