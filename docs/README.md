# TAMA (Rust) documentation

Start with the [**Overview**](OVERVIEW.md) — what TAMA is for, where it fits in an
Iso-Seq pipeline, and a one-paragraph description of every tool.

Then dive into the per-area reference pages:

| Page | Covers |
|---|---|
| [File formats](file-formats.md) | The TAMA **BED12** dialect, the collapse report files, filelists |
| [Collapse](collapse.md) | `tama collapse` — reads → transcript models |
| [Merge](merge.md) | `tama merge` — combine annotations across sources |
| [Filtering](filtering.md) | `tama filter single-read / fragments / polya / primary-orf` |
| [ORF & NMD](orf-nmd.md) | `tama orf seek / blastp-parse / add-cds / extract-cds`, `format bed2gtf-orf` |
| [Read support](read-support.md) | `tama support levels / collapse-cluster / merge-collapse` |
| [File statistics](stats.md) | `tama stats degradation / saturation / model-changes` |
| [Format conversion](format-conversion.md) | `tama format bed2gtf / gtf2bed / gff2bed / fastq2fasta / id-filter` |
| [Utilities](utilities.md) | `tama cleanup polya`, `tama split fasta / sam`, `tama variants call` |
| [Benchmarks](benchmarks.md) | `collapse` speed & memory vs. the original Python 2 |

Each page describes purpose, inputs, outputs (exact filenames), and the important
flags with their defaults. Every flag mirrors the original Python tool's short
option; run `tama <group> <tool> --help` for the complete list.

> These pages document *this Rust port*. For the authoritative biology,
> recommended parameter values, and the original design rationale, see the
> upstream [TAMA wiki](https://github.com/GenomeRIK/tama/wiki) and the
> [TAMA paper](https://doi.org/10.1186/s12864-020-07123-7).
