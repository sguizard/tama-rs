# TAMA (Rust)

A Rust rewrite of [TAMA](https://github.com/GenomeRIK/tama) (Transcriptome
Annotation by Modular Algorithms), a toolkit for processing long-read
(Iso-Seq / Nanopore) transcriptome data.

The original project is a collection of Python 2 scripts. This rewrite ships a
single `tama` binary exposing **every** tool as a subcommand, with native Rust
I/O and **no Python / BioPython / samtools runtime dependency**.

**New to TAMA?** [`docs/OVERVIEW.md`](docs/OVERVIEW.md) explains what the program
is for, the typical Iso-Seq workflow, and the purpose of each tool. The
[`docs/`](docs/README.md) directory has per-tool reference pages (inputs,
outputs, every flag). For the authoritative biology and parameters, see the
original [TAMA wiki](https://github.com/GenomeRIK/tama/wiki).

## Installation

`tama` is compiled from source with Rust's package manager, `cargo`. There is
**no Python / BioPython / samtools** runtime dependency — you only need the Rust
toolchain to build it.

### 1. Install Rust (one time)

```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
. "$HOME/.cargo/env"   # load it into the current shell
```

Accept the defaults when prompted. This installs `cargo` into `~/.cargo/bin`.

### 2. Install `tama`

```sh
cargo install --git https://github.com/sguizard/tama-rs
```

This compiles a release build and installs the `tama` binary to `~/.cargo/bin/`.
To upgrade later, re-run the same command with `--force`.

### 3. Put it on your PATH

If `tama` isn't found afterwards, add `~/.cargo/bin` to your PATH (rustup usually
does this for you):

```sh
echo 'export PATH="$HOME/.cargo/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### 4. Verify

```sh
tama --help
```

## Building from source

To work on the code (or keep the source tree around), clone and build; the binary
lands at `target/release/tama`:

```sh
git clone https://github.com/sguizard/tama-rs.git
cd tama-rs
cargo build --release
./target/release/tama --help
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

By default the tools are quiet, printing only a one-line summary when a step
finishes. Add **`-v` / `--verbose`** (before or after the subcommand) to print
periodic progress to stderr — useful on large inputs to confirm a long
`collapse`/`merge`/`variants` run is still going:

```sh
tama -v merge -f filelist.txt -p merged
# [DEBUG tama::cmd::merge] merge: loaded 120000 transcripts from source ensembl
# [DEBUG tama::cmd::merge] merge: scaffold 3/40 (chr3) done — 5123 genes so far
```

Verbosity can also be set with the `RUST_LOG` environment variable (e.g.
`RUST_LOG=warn` to suppress the summary, `RUST_LOG=off` for silence).

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

## Performance

On `tama collapse` (the compute-heavy core), the Rust port runs **~26–60× faster**
and uses **~5–10× less memory** than the original Python 2, with byte-identical
output. The advantage widens with input size. Full table, methodology, and
caveats: [`docs/benchmarks.md`](docs/benchmarks.md).

## Citation

This project is an independent reimplementation. **If you use it in your
research, you must cite the original TAMA paper:**

> Kuo, R.I., Cheng, Y., Zhang, R., Brown, J.W.S., Smith, J., Archibald, A.L. &
> Burt, D.W. Illuminating the dark side of the human transcriptome with long read
> transcript sequencing. *BMC Genomics* **21**, 751 (2020).
> https://doi.org/10.1186/s12864-020-07123-7

```bibtex
@article{kuo2020tama,
  title   = {Illuminating the dark side of the human transcriptome with long read transcript sequencing},
  author  = {Kuo, Richard I. and Cheng, Yuanyuan and Zhang, Runxuan and Brown, John W. S. and Smith, Jacqueline and Archibald, Alan L. and Burt, David W.},
  journal = {BMC Genomics},
  volume  = {21},
  number  = {1},
  pages   = {751},
  year    = {2020},
  doi     = {10.1186/s12864-020-07123-7}
}
```

Optionally, you may also reference this Rust reimplementation by its repository
URL. See [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.

## License

GPL-3.0-or-later, matching upstream TAMA.
