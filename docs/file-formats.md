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

## Merge output files

`tama merge -p <prefix>` writes:

| File | Contents |
|---|---|
| `<prefix>.bed` | merged transcript models |
| `<prefix>_merge.txt` | per-model provenance: which source model(s) merged into it |
| `<prefix>_trans_report.txt` | per-transcript report incl. `all_source_trans` |
| `<prefix>_gene_report.txt` | per-gene report across sources |

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
