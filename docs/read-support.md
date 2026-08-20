# `tama support` — read support tracking

These tools track exactly *which reads* back each transcript model. The
resulting `_read_support.txt` "levels" file is what the [filters](filtering.md)
and [`stats`](stats.md) tools consume to make support-based decisions.

## `support levels`

Build a per-model read-support file by combining the `trans_read.bed` files from
one or more collapse runs, optionally following a merge.

```sh
tama support levels -f filelist.txt -m merged_merge.txt -o merged
```

- `-f` filelist, one source per line, **3 tab-separated columns**:
  `source_name<TAB>trans_read_file<TAB>file_type` where `file_type` is
  `trans_read` (a collapse `_trans_read.bed`) or `ref_anno`.
- `-m` the merge `_merge.txt` from [`tama merge`](merge.md), or `no_merge` if the
  input isn't merged.
- `-mt` (`tama`|`cupcake`|`filter`) — merge-file type (default `tama`).
- `-o` output prefix → writes `<prefix>_read_support.txt`.

The output records, per merged model, how many reads and how many sources support
it, and the read IDs — the raw material for filtering singletons and poly-A
artefacts.

## `support collapse-cluster`

Read support for collapse models when your reads came through a clustering step
(e.g. Iso-Seq cluster), mapping cluster IDs back to their member reads.

```sh
tama support collapse-cluster ...   # see --help
```

Takes a collapse `_trans_read.bed` and a cluster file; the cluster-file type is
selectable (including `no_cluster`).

## `support merge-collapse`

Aggregate per-source collapse support into merged-model support — the bridge
between individual collapse runs and a merged annotation when you're not using
`levels`.

---

All three preserve read/cluster membership faithfully; only the *order* of read
IDs listed within a cell may differ from the original Python (an artefact of
Python-2 dict ordering) — the membership itself is identical.
