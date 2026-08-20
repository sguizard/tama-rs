# `tama filter`

Clean up an annotation using read support and model geometry. Four filters,
usually applied in sequence.

## `filter single-read`

Remove models supported by only a few reads/sources — likely noise.

```sh
tama filter single-read -b annotation.bed -r read_support.txt -o out
```

- `-b` annotation BED, `-r` [read-support levels file](read-support.md),
  `-o` output prefix.
- `-l` (`gene`|`transcript`) — remove at gene or transcript level (default `gene`).
- `-k` (`keep_multi`|`remove_multi`) — whether to spare multi-exon models
  (default `keep_multi`; multi-exon models are more likely real even with one read).
- `-s` — minimum number of supporting **sources** (default 1).
- `-n` — minimum number of supporting **reads** (default 2).

## `filter fragments`

Remove 5'/3'-**degraded fragment models** that are wholly contained within a
longer model of the same gene, extending the retained model's bounds to cover
them.

```sh
tama filter fragments -f annotation.bed -o out
```

- `-f` input BED, `-o` output prefix.
- `-m` — exon/splice-junction wobble threshold (default 10).
- `-e` — transcript-ends wobble threshold (default 500) — fragments can differ
  substantially at their ends and still be absorbed.
- `-s` — single-exon overlap percent threshold (default 20).

## `filter polya`

Remove models whose supporting reads are genomic **poly-A run-on** artefacts
(genome-templated A-tails rather than true 3' ends).

```sh
tama filter polya -b annotation.bed -f filelist.txt -r read_support.txt -o out
```

- `-b` annotation BED, `-r` read-support file, `-o` prefix.
- `-f` filelist mapping `source_name<TAB>polya_file` (the `_polya.txt` from each
  collapse run).
- `-p` — percent poly-A threshold to call a read an artefact (default 75.0).
- `-l` (`gene`|`transcript`) removal level (default `gene`).
- `-a` (`all_polya`|`singleton_polya`) — remove any poly-A-supported model, or
  only those supported *solely* by poly-A reads (default `singleton_polya`).
- `-k` (`keep_multi`|`remove_multi`) multi-exon handling (default `remove_multi`).

Outputs include `<prefix>_polya_report.txt`, `_polya_support.txt`, and a
`_trash_polya.bed` of what was removed.

## `filter primary-orf`

Keep a single **primary transcript per gene**, chosen by ORF quality then length
— useful for building a one-isoform-per-gene reference.

```sh
tama filter primary-orf -b orf_nmd.bed -o out
```

- `-b` an [ORF/NMD-annotated BED](orf-nmd.md), `-o` output prefix.

Ranking prefers full-length, BLASTP-supported, non-NMD models, breaking ties by
transcript length. Outputs the primary BED plus a `_singleton.bed` /
`_singleton_report.txt` for genes that had only one model, and a `_discarded.txt`.

## Typical order

```sh
tama filter single-read -b merged.bed -r support.txt -o f1
tama filter fragments    -f f1.bed -o f2
tama filter polya        -b f2.bed -f polya_filelist.txt -r support.txt -o f3
# (primary-orf later, after ORF annotation)
```
