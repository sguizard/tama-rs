# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] — 2026-08-21

### Added

- **`tama completions <shell>`.** Generates a completion script on stdout for
  `fish`, `bash`, `zsh`, `powershell`, or `elvish` from the clap definitions, so
  it can never drift from the flags the binary accepts:
  `tama completions fish > ~/.config/fish/completions/tama.fish`. A pre-generated
  fish script is also committed at `completions/tama.fish` (regenerate with
  `completions/regenerate.sh`; CI fails if it goes stale).

### Changed

- **Fixed-value flags are validated at parse time.** The options that accept a
  known set of words — `collapse`/`variants call` `-x` and `-icm`,
  `collapse`/`merge` `-e` and `-d`, `collapse` `-sj` and `-rm`, `filter polya`
  and `filter single-read` `-l`/`-k`, `filter polya -a`, `orf extract-cds`/
  `orf add-cds` `-s`, `orf blastp-parse -f`, `support levels -mt`, and
  `format id-filter` `-f`/`-s` — are now typed rather than free-form strings.
  Accepted spellings are unchanged (`no_cap`, `common_ends`, `singleton_polya`,
  …), and they now tab-complete. **Behaviour change:** a value outside the set is
  a usage error (exit 2) listing the valid options. Previously some of these were
  silently ignored and fell through to the default branch.

## [0.2.1] — 2026-08-21

### Added

- **`-v` / `--verbose` global flag.** Off by default (only end-of-run summaries
  are printed, unlike the original which is very chatty). When set, the
  compute-heavy tools (`collapse`, `merge`, `variants`) print periodic progress
  heartbeats to stderr so you can tell a long run is still going. Also
  controllable via `RUST_LOG`.

## [0.2.0] — 2026-08-21

### Added

- **`tama merge`: no_cap and mixed sources.** `merge` now accepts `no_cap`
  sources and any mix of capped/no_cap in one filelist, porting the original's
  phased grouping (`hunter_prey_capped` → `hunter_prey_mixed` →
  `hunter_prey_nocap`) and the `compare_transcripts_capped_nocap` /
  `_both_nocap` comparisons. A 5'-degraded no_cap model attaches to the capped
  model(s) it matches without merging distinct capped models together.
- **`tama merge -s` / `-cds`.** `-s <source[,…]>` appends a source's original
  `gene_id;transcript_id` to the merged model's column-4 ID; `-cds <source[,…]>`
  copies that source's CDS (thick coords). Both previously errored as
  unimplemented.
- Golden tests for no_cap and mixed merge against the original Python 2
  (`tests/golden_merge_nocap/`).
- **Citation metadata**: `CITATION.cff` and a README Citation section requiring
  users to cite the original TAMA paper (Kuo et al., *BMC Genomics* 21:751, 2020).
- **Benchmarks** ([`docs/benchmarks.md`](docs/benchmarks.md)): `collapse` is
  ~26–60× faster than the Python 2 original with ~5–10× less memory, output
  byte-identical.

### Notes

- `.bed`, `_merge.txt`, and `_gene_report.txt` are byte-identical to the original
  for capped, no_cap, and mixed merges. In `_trans_report.txt`, when a coordinate
  is contended by a capped vs a no_cap member the support-cell "winner" (and, for
  `-s`, the order of appended IDs when a model has several members from the
  source) follows Python-2 dict iteration order and is not reproducible — the
  same class of artifact already documented for capped merge.

## [0.1.1] — 2026-08-20

Documentation release. No functional changes to any tool.

### Added

- **Documentation** under [`docs/`](docs/README.md): a program
  [overview](docs/OVERVIEW.md) (purpose, Iso-Seq workflow, per-tool summary),
  per-tool reference pages (inputs, outputs, every flag with its default), and
  field-by-field [file-format](docs/file-formats.md) descriptions of the TAMA
  BED12 dialect and the collapse/merge report files.
- **Installation instructions** in the README for users new to Rust
  (`cargo install --git …`, PATH setup, verification).
- This **CHANGELOG**.

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

[0.2.1]: https://github.com/sguizard/tama-rs/releases/tag/v0.2.1
[0.2.0]: https://github.com/sguizard/tama-rs/releases/tag/v0.2.0
[0.1.1]: https://github.com/sguizard/tama-rs/releases/tag/v0.1.1
[0.1.0]: https://github.com/sguizard/tama-rs/releases/tag/v0.1.0
