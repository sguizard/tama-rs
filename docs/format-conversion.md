# `tama format` — format conversion

Move annotations between the TAMA BED12 dialect and the formats other tools and
genome browsers expect.

## `format bed2gtf`

TAMA BED12 → Ensembl-style GTF (exon features, no CDS).

```sh
tama format bed2gtf annotation.bed annotation.gtf
```

For a CDS-aware conversion of an ORF/NMD BED, use
[`format bed2gtf-orf`](orf-nmd.md#format-bed2gtf-orf) instead.

## `format gtf2bed`

Import a GTF into TAMA BED12. The `--source` flag selects the dialect:

```sh
tama format gtf2bed --source ensembl in.gtf out.bed
```

- `--source ensembl | ncbi | stringtie`

Each source lays out its attributes differently; picking the right one ensures
gene/transcript IDs land in TAMA's column-4 `gene_id;transcript_id` form.

## `format gff2bed`

Import a GFF3 into TAMA BED12.

```sh
tama format gff2bed --source cupcake in.gff out.bed
```

- `--source cupcake | liftoff`

## `format fastq2fasta`

Convert Nanopore FASTQ reads to FASTA (dropping quality) prior to mapping.

```sh
tama format fastq2fasta reads.fastq reads.fasta
```

## `format id-filter`

Restructure or trim the BED column-4 ID field — e.g. promote an Ensembl ID,
reorder the `;`-separated subfields, or drop the ORF flags before handing the BED
to another program.

```sh
tama format id-filter -b in.bed -o out.txt -s ensembl_merge
```

- `-b` input BED, `-o` output.
- `-s` sub-field method: `ensembl_merge`, `ensembl_orf`, or `custom`.
- `-r` custom reshuffle spec, e.g. `3,4,1,2` (with `-s custom`).
- `-d` sub-field delimiters (default `;`).
- `-f` filter level: `none` or `only_match`.

See [File formats](file-formats.md#the-column-4-id-field) for what the subfields
mean.
