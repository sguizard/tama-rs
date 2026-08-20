# Utilities: cleanup, split, variants

## `tama cleanup polya`

Trim genomic **poly-A tails** off FLNC reads *before* mapping, so that
genome-templated A-stretches don't get mistaken for real transcript 3' ends
downstream.

```sh
tama cleanup polya -f flnc.fa -p cleaned -m 200
```

- `-f` FLNC FASTA, `-p` output prefix.
- `-m` minimum read length to keep after trimming (default 200).

Outputs the cleaned reads plus reports:
`<prefix>_polya_flnc_report.txt`, `_tails.fa` (the removed tails),
`_discarded_reads.txt` (too short after trimming), and `_summary.txt`.

## `tama split`

Split large inputs so heavy steps can run in parallel.

### `split fasta`

```sh
tama split fasta transcripts.fa chunk 8
```

Positional args: `<fasta> <prefix> <split_number>`. Writes `chunk_1.fa` …
`chunk_8.fa` — e.g. to fan a big protein FASTA out across parallel BLASTP jobs.

### `split sam`

```sh
tama split sam mapped.sam 4 part
```

Positional args: `<sam> <num_files> <prefix>`. Splits a mapped SAM by chromosome
into N files (`part_1.sam` …), for parallel collapse.

## `tama variants call`

Call sequence variants (mismatches, indels, soft-clips) that are supported by
enough reads at a genomic position, with per-position read coverage. It reuses
[`collapse`](collapse.md)'s per-read CIGAR/variation analysis — the variant
outputs depend only on per-read variation and coverage, not on the collapse
grouping, so this is a standalone entry point for the same logic.

```sh
tama variants call -s sorted.sam -f genome.fa -p variants -x no_cap
```

- `-s` sorted SAM (or BAM with `-b`), `-f` genome FASTA, `-p` output prefix.
- `-x` `capped`|`no_cap`, `-c` coverage %, `-i` identity %, `-icm` identity method
  — same meanings as [`collapse`](collapse.md#key-parameters).
- `-sjt` splice-junction error threshold (default 10).
- `-vc` variation coverage threshold — the minimum number of reads that must show
  a variant for it to be reported (default 5).

Outputs `<prefix>_variants.txt` (the called variants) and `<prefix>_varcov.txt`
(read coverage per variant position).

> Only the order of read IDs inside the `cluster_list` column may differ from the
> original Python (a dict-ordering artefact); the called variants and coverage
> are identical.
