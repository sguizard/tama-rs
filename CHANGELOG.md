# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] — 2026-08-20

Initial release: a complete Rust rewrite of
[TAMA](https://github.com/GenomeRIK/tama) (Transcriptome Annotation by Modular
Algorithms), the long-read (Iso-Seq / Nanopore) transcriptome annotation toolkit.

### Added

- **All 34 original TAMA tools**, exposed as subcommands of a single `tama`
  binary:
  - Core: `collapse`, `merge`.
  - `filter` — `single-read`, `fragments`, `polya`, `primary-orf`.
  - `orf` — `seek`, `blastp-parse`, `add-cds`, `extract-cds`.
  - `support` — `levels`, `collapse-cluster`, `merge-collapse`.
  - `stats` — `degradation`, `saturation`, `model-changes`.
  - `format` — `bed2gtf`, `bed2gtf-orf`, `gtf2bed` (ensembl/ncbi/stringtie),
    `gff2bed` (cupcake/liftoff), `fastq2fasta`, `id-filter`.
  - `cleanup polya`, `split` (`fasta`, `sam`), `variants call`.
- **Native Rust I/O** — no Python, BioPython, or samtools runtime dependency;
  the toolkit builds to one self-contained binary.
- **Documentation** under [`docs/`](docs/README.md): a program
  [overview](docs/OVERVIEW.md), per-tool reference pages, and field-by-field
  [file-format](docs/file-formats.md) descriptions.
- **Installation instructions** in the README for users new to Rust
  (`cargo install --git …`).
- Cargo workspace: `tama-core` (model + algorithms), `tama-io` (format I/O),
  `tama` (CLI).
- CI (`.github/workflows/ci.yml`): rustfmt, clippy (`-D warnings`), build, test.
- GPL-3.0-or-later license, matching upstream TAMA.

### Compatibility

Each tool is validated against the original Python 2 script by running both on
the same input and diffing (“golden” tests). Outputs are **byte-identical** to
the original, with one class of exception: a few *diagnostic report columns*
list the same tokens in a different order, because the original derives that
order from Python-2 dict iteration (hash order), which is neither reproducible
nor semantically meaningful. These cells are compared as sets in the tests.
Affected: `collapse` `_trans_report` `collapse_error_nuc`; `merge`
`all_source_trans`; `support` read/cluster lists; `variants` `cluster_list`;
`stats model-changes` diff cells.

### Known limitations

- **`collapse`**: BAM input (`-b`), multimap handling, the
  `_local_density_error.txt` output, and the `-rm low_mem` run mode are
  recognised but not yet fully ported. SAM input in the default run mode is the
  validated path. Per-model variant output (`_variants.txt` / `_varcov.txt`) is
  produced by `tama variants call`, which reuses the same per-read analysis.
- **`merge`**: only `capped` sources are supported so far; `no_cap`/mixed-source
  merging and the `-s`/`-cds` source overrides are not yet ported.

[0.1.0]: https://github.com/sguizard/tama-rs/releases/tag/v0.1.0
