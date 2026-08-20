# `tama collapse`

Build non-redundant transcript models from mapped long reads.

```sh
tama collapse -s sorted.sam -f genome.fa -p out -x no_cap
```

## What it does

`collapse` reads a coordinate-sorted alignment file and, for each read:

1. Walks the **CIGAR** to reconstruct exon boundaries and measure mapping
   coverage and identity.
2. Rejects reads below the coverage/identity cutoffs.
3. Detects **genomic poly-A run-on** (the read's 3' end sits on a genome-encoded
   A-stretch, so the "transcript end" is an artefact) and **strand conflicts**
   (SAM-flag strand vs. the `XS:A:` splice-site tag).
4. Groups surviving reads into genes by locus overlap, then **collapses** reads
   that represent the same isoform into one model, choosing representative exon
   boundaries within the wobble thresholds.

The result is one transcript model per distinct isoform, each with a full record
of its supporting reads.

## Capped vs. no_cap (`-x`)

This is the most important choice.

- **`capped`** — for 5'-cap-selected libraries, where the 5' end is trustworthy.
  Two reads only collapse together if they have the **same exon count**; a
  shorter read is treated as a genuinely different (shorter) transcript.
- **`no_cap`** (default) — tolerates 5' RNA **degradation**. A read with fewer
  exons or a shorter first exon is accepted as *support for a longer model*, as
  long as its 3' end and internal splice junctions match. Use this when 5' ends
  are unreliable.

Running both and comparing them is exactly what [`stats degradation`](stats.md)
does to quantify a library's degradation level.

## Inputs

- `-s, --sam` — coordinate-sorted SAM (or BAM with `-b`). Read as plain TSV:
  query name, flag, reference, position, CIGAR, SEQ, and the optional `XS:A:`
  strand tag.
- `-f, --fasta` — genome FASTA (used for poly-A run-on detection and variation).
- `-p, --prefix` — output prefix.

## Outputs

`<prefix>.bed`, `_trans_read.bed`, `_read.txt`, `_trans_report.txt`,
`_polya.txt`, `_strand_check.txt`, `_report.txt` — see
[File formats](file-formats.md#collapse-output-files).

## Key parameters

| Flag | Default | Meaning |
|---|---|---|
| `-x` | `no_cap` | `capped` or `no_cap` (above) |
| `-e` | `common_ends` | exon-end selection: `common_ends` or `longest_ends` |
| `-c` | `99.0` | minimum mapping **coverage** percent to accept a read |
| `-i` | `85.0` | minimum **identity** percent to accept a read |
| `-icm` | `ident_cov` | identity calc method: `ident_cov` or `ident_map` |
| `-a` | `10` | **5'** wobble threshold (bp) for calling two starts "the same" |
| `-m` | `10` | **exon / splice-junction** wobble threshold (bp) |
| `-z` | `10` | **3'** wobble threshold (bp) |
| `-d` | `merge_dup` | duplicate handling: `merge_dup` or `no_merge` |
| `-sj` | `no_priority` | splice-junction priority: `no_priority` or `sj_priority` |
| `-sjt` | `10` | splice-junction error threshold (bp) |
| `-lde` | `1000` | local density error threshold |
| `-vc` | `5` | variation coverage threshold (reads) |
| `-rm` | `original` | run mode: `original` or `low_mem` |
| `-b` | off | treat input as BAM |

`-e longest_ends` keeps the longest observed 5'/3' end for a model instead of the
most common one. The `-a`/`-m`/`-z` thresholds are the tolerance for treating two
coordinates as identical — larger values collapse more aggressively.

## Notes and current limitations

- BAM input (`-b`), multimap handling, the `_local_density_error.txt` output, and
  `-rm low_mem` are recognised but not yet fully ported — SAM input with the
  default run mode is the validated path.
- Per-model variant/variation output (`_variants.txt`/`_varcov.txt`) is produced
  by [`tama variants call`](utilities.md#tama-variants-call), which reuses the same
  per-read analysis.
