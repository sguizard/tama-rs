# File formats

TAMA tools pass data around as a small number of text formats. This page
describes the ones you'll touch most: the TAMA **BED12** dialect and the report
files that `collapse`/`merge` emit.

## TAMA BED12

A standard 12-column BED, with TAMA-specific conventions in a few columns.

| Col | BED name | TAMA usage |
|----:|----------|------------|
| 1 | chrom | scaffold / chromosome |
| 2 | chromStart | model start (0-based) |
| 3 | chromEnd | model end |
| 4 | name | **`gene_id;transcript_id`** (see below) |
| 5 | score | unused (often `40`) |
| 6 | strand | `+` / `-` |
| 7 | thickStart | CDS start (= chromStart when no CDS) |
| 8 | thickEnd | CDS end (= chromEnd when no CDS) |
| 9 | itemRgb | source colour (merge) or `255,0,0` |
| 10 | blockCount | number of exons |
| 11 | blockSizes | comma-separated exon lengths |
| 12 | blockStarts | comma-separated exon offsets from chromStart |

### The column-4 ID field

The `name` column packs several `;`-separated subfields. At minimum:

```
G1;G1.1
└┬┘ └─┬─┘
gene  transcript
```

Downstream tools **append** subfields to column 4. After ORF/NMD annotation
(`orf add-cds`) it becomes:

```
gene_id;transcript_id;prot_id;degrade_flag;match_flag;nmd_flag;frame
```

- `degrade_flag` — `full_length` vs `prot_ok`/degraded (5' truncation detected).
- `match_flag` — whether the ORF had BLASTP support.
- `nmd_flag` — NMD candidate (premature stop upstream of the last splice junction).
- `frame` — reading frame used.

The `format id-filter` tool exists precisely to reshuffle/trim these subfields
when you need to hand the BED to another program.

## Collapse output files

`tama collapse -p <prefix>` writes:

| File | Contents |
|---|---|
| `<prefix>.bed` | the collapsed transcript models (TAMA BED12) |
| `<prefix>_trans_read.bed` | every input read, in BED12, tagged with the model it supports |
| `<prefix>_read.txt` | per-read report: mapping flags, coverage/identity, strand, classification |
| `<prefix>_trans_report.txt` | per-model report: supporting read count, error stats, source reads |
| `<prefix>_polya.txt` | reads flagged as genomic poly-A run-on |
| `<prefix>_strand_check.txt` | reads whose SAM flag strand disagrees with the `XS:A:` tag |
| `<prefix>_report.txt` | run summary (counts by category) |

> Not yet emitted by `collapse` itself: `_variants.txt` / `_varcov.txt` (use
> `tama variants call`) and `_local_density_error.txt`.

### `_read.txt` columns

Per input read, one line (header row included):

| Column | Meaning |
|---|---|
| `read_id` | query name from the SAM |
| `mapped_flag` | strand/mapping category (`forward_strand`, `reverse_strand`, `unmapped`, `not_primary`, `chimeric`, `unknown`) |
| `accept_flag` | whether the read passed the coverage/identity cutoffs |
| `percent_coverage` | fraction of the read aligned |
| `percent_identity` | aligned-base identity |
| `error_line<h;s;i;d;m>` | error tallies: hard-clip; soft-clip; insertion; deletion; mismatch |
| `length` | mapped length |
| `cigar` | the read's CIGAR |

### `_trans_report.txt` columns

Per collapsed model:

| Column | Meaning |
|---|---|
| `transcript_id` | model ID |
| `num_clusters` | number of supporting reads/clusters |
| `high_coverage` / `low_coverage` | supporting-read counts above/below the coverage cutoff |
| `high_quality_percent` / `low_quality_percent` | corresponding percentages |
| `start_wobble_list` / `end_wobble_list` | spread of observed 5'/3' ends across supporting reads |
| `collapse_sj_start_err` / `collapse_sj_end_err` | splice-junction error at the chosen boundaries |
| `collapse_error_nuc` | per-nucleotide error detail (token order within this cell is not significant — see the equivalence note in the README) |

### `_polya.txt` columns

Reads flagged as genomic poly-A run-on:

`cluster_id`, `trans_id`, `strand`, `a_percent` (fraction of the tail window that
is A), `a_count`, `sequence`.

### `_strand_check.txt` columns

Reads whose SAM-flag strand disagrees with the `XS:A:` splice tag:
`read_id`, `scaff_name`, `start_pos`, `cigar`, `strands`.

### `_report.txt`

Key/value run summary: `Total Gene Count`, `Total Transcript Count`,
`Total Accepted Reads`, `Total Discarded Reads`.

### `_trans_read.bed`

Standard TAMA BED12 — one line **per input read**, with column 4 set to the
`gene_id;transcript_id` of the model that read supports. This is the file
`support levels` and `stats degradation` consume.

## Merge output files

`tama merge -p <prefix>` writes:

| File | Contents |
|---|---|
| `<prefix>.bed` | merged transcript models |
| `<prefix>_merge.txt` | per-model provenance: which source model(s) merged into it |
| `<prefix>_trans_report.txt` | per-transcript report incl. `all_source_trans` |
| `<prefix>_gene_report.txt` | per-gene report across sources |

### merge `_trans_report.txt` columns

`transcript_id`, `num_clusters`, `sources`, `start_wobble_list`,
`end_wobble_list`, `exon_start_support`, `exon_end_support`, `all_source_trans`
(the list of contributing per-source transcript IDs; token order within the cell
is not significant).

### merge `_gene_report.txt` columns

`gene_id`, `num_clusters`, `num_final_trans`, `sources`, `chrom`, `start`, `end`,
`source_genes`, `source_summary`.

### merge `_merge.txt`

No header — TAMA BED12 lines for each **source** model, tagged with the merged
model it folded into, so you can trace provenance model-by-model.

## Read-support levels file

`support levels` writes `<prefix>_read_support.txt` with one line per merged model:

| Column | Meaning |
|---|---|
| `merge_gene_id` / `merge_trans_id` | merged model IDs |
| `gene_read_support` / `trans_read_support` | supporting-read counts at gene / transcript level |
| `source_prefix` | contributing source name(s) |
| `source_trans_line` | contributing per-source transcript IDs |
| `source_read_line` | contributing read IDs |

This is the file the [filters](filtering.md) and
[`stats saturation` / `model-changes`](stats.md) read.

## Filelists

Two tools take a tab-separated *filelist* describing multiple sources.

**`tama merge -f`** — 4 columns per line:

```
<bed_file>	<capped|no_cap>	<start,junction,end priority>	<source_id>
```

e.g. `sampleA.bed	capped	1,1,1	sampleA`. Lower priority numbers win when
choosing coordinates; the three numbers set start / splice-junction / end
priority independently.

**`tama support levels -f`** — 3 columns per line:

```
<source_name>	<trans_read.bed file>	<file_type>
```

where `file_type` is `trans_read` (collapse output) or `ref_anno`.

See each tool's page for how these are consumed.
