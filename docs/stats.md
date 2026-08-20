# `tama stats` — file statistics

Diagnostic reports about a library or an annotation.

## `stats degradation`

Quantify how much 5' RNA **degradation** is in a library by comparing a
**capped** collapse run against a **no_cap** run of the same data (the
"degradation signature"). Lots of models that only appear in `no_cap` — or
`no_cap` models with many more supporting reads — indicate heavy 5' degradation.

```sh
tama stats degradation -c capped_trans_read.bed -nc nocap_trans_read.bed -o degradation.txt
```

- `-c` the **capped** collapse `trans_read.bed`.
- `-nc` the **no_cap** collapse `trans_read.bed`.
- `-o` output report (includes the degradation signature value and per-category
  transcript/read counts).

See [`collapse -x`](collapse.md#capped-vs-no_cap--x) for what capped/no_cap mean.

## `stats saturation`

A sampling **saturation curve**: how many genes/transcripts are discovered as a
function of reads sequenced. A curve that has plateaued means more sequencing
won't find much new; one still climbing means you're undersampled.

```sh
tama stats saturation -r read_support.txt -b 1000 -o saturation.txt
```

- `-r` a [read-support levels file](read-support.md).
- `-b` read bin size (points along the curve).
- `-o` output table.

## `stats model-changes`

Find reads/transcripts that map to **different genes or transcripts** between two
sources — e.g. a reference annotation vs. an alternative one. Useful for spotting
where two annotations disagree about gene structure.

```sh
tama stats model-changes -b annotation.bed -r read_support.txt --ref refname --alt altname -o out
```

- `-b` annotation BED, `-r` read-support file, `-o` output prefix.
- `-ref` / `-alt` — the two source names to compare (default `NA`).

Writes several diff reports: `<prefix>_diff_report.txt`, `_diff_genes.txt`,
`_diff_trans.txt`, and the `_diff_one_source_*` variants for models present in
only one source.

> A meaningful run needs a **multi-source** read-support file; comparing a single
> source produces only headers.
