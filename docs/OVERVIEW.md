# TAMA — What it does and how the tools fit together

TAMA (**T**ranscriptome **A**nnotation by **M**odular **A**lgorithms) builds
transcriptome/genome annotations from **long-read RNA sequencing** (PacBio
Iso-Seq, Oxford Nanopore), and works with other data types too. This document
explains the *purpose* of each tool and how they chain into a workflow. It
describes the concepts, not every flag — run `tama <group> <tool> --help` for
the exact parameters, and see the original
[TAMA wiki](https://github.com/GenomeRIK/tama/wiki) and
[paper](https://doi.org/10.1186/s12864-020-07123-7) for the authoritative
biology and recommended settings.

## The problem TAMA solves

When you sequence full-length transcripts and map them to a genome, you get
**many reads per gene that are slight variants of each other** — the same
biological isoform showing up repeatedly with small differences from sequencing
error, RNA degradation (5' truncation), or imprecise transcript ends. Left as-is
this inflates transcript counts and muddies every downstream analysis.

TAMA's job is to turn that pile of mapped reads into a **clean, non-redundant set
of transcript models** while keeping a full record of *which reads support which
model* — so you can trust the annotation and trace it back to the evidence.

## Where TAMA sits in an Iso-Seq pipeline

```
subreads → CCS → FLNC (adapter/poly-A/concatemer removal) → [cluster/polish]
        → map to genome (minimap2 / GMAP)  →  sorted SAM/BAM
        → TAMA COLLAPSE           (per-sample transcript models)
        → TAMA MERGE              (combine samples/methods, keep provenance)
        → TAMA GO tools           (filter, annotate ORFs/NMD, read support, stats, format)
```

The two **core** tools are `collapse` and `merge`; the **TAMA GO** tools are a
suite of smaller, composable utilities you apply afterwards.

## The interchange format: TAMA BED12

Almost every tool reads/writes a **BED12** file with one TAMA-specific
convention: column 4 (the `name` field) packs identifiers as
`gene_id;transcript_id` (and some tools append more subfields, e.g. ORF/NMD adds
`;prot_id;degrade_flag;match_flag;nmd_flag;frame`). Columns 7/8
(`thickStart`/`thickEnd`) carry the CDS region when present. This shared format
is why the tools compose so freely.

---

## Core tools

### `tama collapse` — build transcript models from mapped reads
Groups mapped reads by genomic locus and **collapses redundant reads into single
transcript models**, choosing representative exon boundaries. It also flags
low-quality mappings, detects **genomic poly-A run-on** (reads that ran into a
genome-encoded A-stretch rather than a real transcript end), checks for
**strand ambiguity**, and records per-read support.

- **Input:** sorted SAM/BAM + genome FASTA. **Output:** `<prefix>.bed` plus
  read/transcript reports, poly-A, strand-check, and (with the variant tool)
  variation files.
- **`-x capped` vs `-x no_cap`:** *capped* libraries used 5'-cap selection, so a
  read's 5' end is trustworthy — models must share the same exon count to merge.
  *no_cap* accepts 5'-**degraded** reads (fewer exons / shorter first exon) as
  support for a longer model, as long as the 3' end and internal splice
  junctions agree.
- **Key thresholds:** minimum coverage (99%) and identity (85%) to accept a read;
  5'/exon/3' wobble tolerances (10 bp) for calling two ends "the same".

### `tama merge` — combine annotations across sources
Takes several BED annotations (e.g. different samples, or Iso-Seq + a reference)
and produces **one merged annotation that records which source each model came
from**. Per-source **priorities** let you decide whose start/junction/end
coordinates win when models are merged. Merge uses looser 5'/3' tolerances
(20 bp) than collapse by default.

---

## TAMA GO — the utility suite

### Transcript filtering (`tama filter …`)
Clean up an annotation using read support and model geometry:
- **`single-read`** — remove models supported by only one read (likely noise),
  at gene or transcript level, with options to keep multi-exon models.
- **`fragments`** — remove 5'/3'-**degraded fragment models** that are contained
  within a longer model of the same gene (extends the retained model's bounds).
- **`polya`** — remove models whose reads are poly-A run-on artefacts
  (genome-templated A-tails rather than true 3' ends).
- **`primary-orf`** — keep one **primary transcript per gene**, ranked by ORF
  quality (full-length, blastp match, no NMD), then length.

### ORF / NMD prediction (`tama orf …`, `tama format bed2gtf-orf`)
Predict coding regions and nonsense-mediated-decay status:
- **`seek`** — find open reading frames in transcript sequences (all 3 frames),
  emitting candidate proteins for homology search.
- **`blastp-parse`** — parse BLASTP results to pick the best-supported ORF per
  transcript.
- **`add-cds`** — write the chosen CDS coordinates back onto the BED and flag
  **NMD** candidates (stop codon well upstream of the last splice junction) and
  **5'-degraded** models.
- **`extract-cds`** — output just the CDS region of each coding model.
- **`format bed2gtf-orf`** — export an ORF/NMD BED to GTF with
  exon/CDS/UTR/start/stop features.

### Read support (`tama support …`)
Track exactly which reads back each model — essential for allele-specific and
quantitative work:
- **`collapse-cluster`** — read support for collapse models from a clustering
  file.
- **`levels`** — build a per-model read-support file from pre-merge
  `trans_read.bed` files (+ an optional merge file).
- **`merge-collapse`** — aggregate per-source collapse support into
  merged-model support.

### File statistics (`tama stats …`)
- **`degradation`** — compare a **capped** vs **no_cap** run to quantify how much
  5' degradation is in your library (the "degradation signature").
- **`saturation`** — a sampling curve of genes discovered vs. reads sequenced
  (how close to saturation your sequencing depth is).
- **`model-changes`** — find reads that map to **different genes/transcripts**
  between two sources (e.g. reference vs. alternative annotation).

### Sequence cleanup (`tama cleanup polya`)
Trim genomic **poly-A tails** off FLNC reads before mapping, classifying reads as
passing / too-short / over-trimmed and reporting tail statistics.

### Format conversion (`tama format …`)
Interoperate with other tools and browsers:
- **`bed2gtf`** / **`bed2gtf-orf`** — TAMA BED → Ensembl-style GTF (with CDS for
  the ORF variant).
- **`gtf2bed --source ensembl|ncbi|stringtie`**,
  **`gff2bed --source cupcake|liftoff`** — import annotations from other tools
  into TAMA BED12.
- **`fastq2fasta`** — Nanopore FASTQ → FASTA.
- **`id-filter`** — restructure/filter the BED ID field (e.g. promote Ensembl IDs).

### Splitting (`tama split …`)
Parallelize heavy steps:
- **`fasta`** — split a FASTA into N chunks (e.g. for parallel BLASTP).
- **`sam`** — split a mapped SAM by chromosome into N files.

### Variant calling (`tama variants call`)
Call sequence variants (mismatches, indels, clips) supported by enough reads at a
position, with per-position read coverage — reusing collapse's per-read analysis.
Outputs `_variants.txt` and `_varcov.txt`.

---

## A typical end-to-end run

```sh
# 1. Per-sample models
tama collapse -s sampleA.sorted.sam -f genome.fa -p sampleA -x no_cap

# 2. Merge samples (filelist.txt lists each sample's BED, cap flag, priority, name)
tama merge -f filelist.txt -p merged

# 3. Read support, then filter out singletons and fragments
tama support levels -f trans_read_filelist.txt -m merged_merge.txt -o merged
tama filter single-read -b merged.bed -r merged_read_support.txt -o merged_f1
tama filter fragments   -f merged_f1.bed -o merged_f2

# 4. ORF / NMD annotation
tama orf seek -f merged_f2_transcripts.fa -o orfs.faa
#   ... run blastp on orfs.faa ...
tama orf blastp-parse -b blastp.out -o blastp_parse.txt
tama orf add-cds -p blastp_parse.txt -a merged_f2.bed -f merged_f2_transcripts.fa -o merged_orf.bed
tama format bed2gtf-orf merged_orf.bed merged_orf.gtf
```

For the biological rationale, recommended parameters, and file-format details,
consult the original [TAMA wiki](https://github.com/GenomeRIK/tama/wiki).
