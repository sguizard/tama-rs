# `tama orf` — ORF & NMD prediction

Predict coding regions and nonsense-mediated-decay (NMD) status for transcript
models. The workflow is: find ORFs → search them against a protein DB with
BLASTP (external) → parse the hits → write CDS coordinates back onto the BED.

```
transcripts.fa ── orf seek ──► orfs.faa ── (blastp) ──► blastp.out
    orf blastp-parse ──► blastp_parse.txt
    orf add-cds (+ annotation.bed, transcripts.fa) ──► orf_nmd.bed
    format bed2gtf-orf ──► orf_nmd.gtf
```

## `orf seek`

Find open reading frames in transcript sequences (all three frames) and emit
candidate proteins for homology search.

```sh
tama orf seek -f transcripts.fa -o orfs.faa
```

- `-f` transcript FASTA, `-o` output protein FASTA.

## `orf blastp-parse`

Parse a BLASTP **pairwise text** output and pick the best-supported ORF per
transcript.

```sh
tama orf blastp-parse -b blastp.out -o blastp_parse.txt
```

- `-b` BLASTP pairwise output, `-o` parse file.
- `-f` (`uniref`|`ensembl`) — how to interpret the DB subject IDs (default `uniref`).

## `orf add-cds`

Write the chosen CDS coordinates back onto the annotation BED, and flag NMD and
5'-degraded models.

```sh
tama orf add-cds -p blastp_parse.txt -a annotation.bed -f transcripts.fa -o orf_nmd.bed
```

- `-p` blastp-parse file, `-a` annotation BED, `-f` transcript FASTA, `-o` output BED.
- `-s` — `include_stop` to keep the stop codon inside the CDS.
- `-d` — distance (bp) from the last splice junction beyond which a premature stop
  triggers the **NMD** flag (default 50).

The output BED's column 4 gains the ORF subfields (`prot_id;degrade_flag;
match_flag;nmd_flag;frame` — see [File formats](file-formats.md#the-column-4-id-field)),
and columns 7/8 carry the CDS thick coordinates.

## `orf extract-cds`

Output just the CDS region of each coding model.

```sh
tama orf extract-cds -b orf_nmd.bed -s no_stop_codon -o cds.bed
```

- `-b` ORF/NMD BED, `-o` output BED.
- `-s` — `include_stop` to keep the stop codon (default excludes it).

## `format bed2gtf-orf`

Export an ORF/NMD BED to GTF with `exon`, `CDS`, `five_prime_utr`,
`three_prime_utr`, `start_codon`, and `stop_codon` features.

```sh
tama format bed2gtf-orf orf_nmd.bed orf_nmd.gtf
```

(This is the CDS-aware sibling of [`format bed2gtf`](format-conversion.md).)
