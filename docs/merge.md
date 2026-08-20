# `tama merge`

Combine several transcript annotations into one, keeping a record of which
source each model came from.

```sh
tama merge -f filelist.txt -p merged
```

## What it does

`merge` takes multiple BED annotations — different samples, different mappers, or
Iso-Seq plus a reference — and produces a **single non-redundant annotation**.
Models that agree (within the end/junction thresholds) are merged into one, and
the `_merge.txt` file records every source model that contributed. Where sources
disagree about a start, splice junction, or 3' end, per-source **priorities**
decide whose coordinate wins.

It reuses the same compare/collapse machinery as [`collapse`](collapse.md), but
across annotations rather than reads.

## The filelist (`-f`)

A tab-separated file, one source per line, **4 columns**:

```
<bed_file>	<capped|no_cap>	<start,junction,end>	<source_id>
```

- **bed_file** — a TAMA BED12 annotation (e.g. a collapse output).
- **capped / no_cap** — how that source's 5' ends should be treated.
- **priority** — three integers for start / splice-junction / 3'-end priority;
  **lower wins**. This lets a trusted reference dictate junctions while an
  Iso-Seq source dictates ends, for instance.
- **source_id** — a short name; used in reports and to colour the output BED.

Example:

```
ref.bed	capped	1,1,1	ensembl
isoseq.bed	no_cap	2,2,2	isoseq
```

## Outputs

| File | Contents |
|---|---|
| `<prefix>.bed` | merged models (source colour in the RGB column) |
| `<prefix>_merge.txt` | which source model(s) merged into each output model |
| `<prefix>_trans_report.txt` | per-transcript report incl. `all_source_trans` |
| `<prefix>_gene_report.txt` | per-gene summary across sources |

## Key parameters

| Flag | Default | Meaning |
|---|---|---|
| `-e` | `common_ends` | exon-end selection: `common_ends` or `longest_ends` |
| `-a` | `20` | **5'** threshold (bp) — note the merge default is looser than collapse's 10 |
| `-m` | `10` | exon / splice-junction threshold (bp) |
| `-z` | `20` | **3'** threshold (bp) — again looser than collapse |
| `-d` | `no_merge` | duplicate handling: `no_merge` or `merge_dup` |
| `-s` | — | take gene/transcript IDs from this source (not yet implemented) |
| `-cds` | — | take CDS from this source (not yet implemented) |

The looser default 5'/3' thresholds (20 vs. collapse's 10) reflect that
independently-built annotations disagree about exact ends more than reads within
one sample do.

## Current limitations

- Only **capped** sources are supported so far; `no_cap`/mixed-source merging
  (`hunter_prey_nocap`/`_mixed` in the original) and the `-s`/`-cds` source
  overrides are not yet ported.
