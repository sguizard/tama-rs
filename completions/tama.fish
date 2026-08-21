# Print an optspec for argparse to handle cmd's options that are independent of any subcommand.
function __fish_tama_global_optspecs
    string join \n v/verbose h/help V/version
end

function __fish_tama_needs_command
    # Figure out if the current invocation already has a command.
    set -l cmd (commandline -opc)
    set -e cmd[1]
    argparse -s (__fish_tama_global_optspecs) -- $cmd 2>/dev/null
    or return
    if set -q argv[1]
        # Also print the command, so this can be used to figure out what it is.
        echo $argv[1]
        return 1
    end
    return 0
end

function __fish_tama_using_subcommand
    set -l cmd (__fish_tama_needs_command)
    test -z "$cmd"
    and return 1
    contains -- $cmd[1] $argv
end

complete -c tama -n "__fish_tama_needs_command" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_needs_command" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c tama -n "__fish_tama_needs_command" -s V -l version -d 'Print version'
complete -c tama -n "__fish_tama_needs_command" -f -a "collapse" -d 'Collapse mapped reads into transcript models (tama_collapse)'
complete -c tama -n "__fish_tama_needs_command" -f -a "merge" -d 'Merge transcript annotations across sources (tama_merge)'
complete -c tama -n "__fish_tama_needs_command" -f -a "format" -d 'Format converters (BED <-> GTF/GFF, FASTQ -> FASTA, ...)'
complete -c tama -n "__fish_tama_needs_command" -f -a "filter" -d 'Filter transcript models'
complete -c tama -n "__fish_tama_needs_command" -f -a "orf" -d 'ORF / NMD prediction tools'
complete -c tama -n "__fish_tama_needs_command" -f -a "support" -d 'Read-support tracking tools'
complete -c tama -n "__fish_tama_needs_command" -f -a "stats" -d 'File statistics tools'
complete -c tama -n "__fish_tama_needs_command" -f -a "split" -d 'File splitting tools'
complete -c tama -n "__fish_tama_needs_command" -f -a "variants" -d 'Variant calling'
complete -c tama -n "__fish_tama_needs_command" -f -a "cleanup" -d 'Sequence cleanup tools'
complete -c tama -n "__fish_tama_needs_command" -f -a "completions" -d 'Generate a shell completion script on stdout'
complete -c tama -n "__fish_tama_needs_command" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand collapse" -s s -l sam -d 'Sorted SAM file (or BAM with --bam). (`-s`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand collapse" -s f -l fasta -d 'Genome FASTA file. (`-f`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand collapse" -s p -l prefix -d 'Output prefix. (`-p`)' -r
complete -c tama -n "__fish_tama_using_subcommand collapse" -s x -l cap-flag -d 'Capped flag. (`-x`)' -r -f -a "capped\t'5\' cap-trimmed data; 5\' ends are trusted'
no_cap\t'Non-cap-trimmed data; 5\' degradation is expected'"
complete -c tama -n "__fish_tama_using_subcommand collapse" -s e -l ends -d 'Collapse exon ends. (`-e`)' -r -f -a "common_ends\t'Use the most common start/end among the collapsed models'
longest_ends\t'Use the longest start/end among the collapsed models'"
complete -c tama -n "__fish_tama_using_subcommand collapse" -s c -l coverage -d 'Minimum coverage percent. (`-c`)' -r
complete -c tama -n "__fish_tama_using_subcommand collapse" -s i -l identity -d 'Minimum identity percent. (`-i`)' -r
complete -c tama -n "__fish_tama_using_subcommand collapse" -l icm -d 'Identity calculation method. (`-icm`)' -r -f -a "ident_cov\t'Include hard/soft clipping in the denominator (default)'
ident_map\t'Exclude hard/soft clipping'"
complete -c tama -n "__fish_tama_using_subcommand collapse" -s a -l five-prime -d '5\' threshold. (`-a`)' -r
complete -c tama -n "__fish_tama_using_subcommand collapse" -s m -l exon-thresh -d 'Exon/splice-junction threshold. (`-m`)' -r
complete -c tama -n "__fish_tama_using_subcommand collapse" -s z -l three-prime -d '3\' threshold. (`-z`)' -r
complete -c tama -n "__fish_tama_using_subcommand collapse" -s d -l dup -d 'Duplicate merge behaviour. (`-d`)' -r -f -a "merge_dup\t'Merge duplicate models into one'
no_merge\t'Keep duplicate models separate'"
complete -c tama -n "__fish_tama_using_subcommand collapse" -l sj -d 'Splice-junction priority. (`-sj`)' -r -f -a "no_priority\t'Do not prioritise splice junctions over ends'
sj_priority\t'Prioritise splice junctions over ends'"
complete -c tama -n "__fish_tama_using_subcommand collapse" -l sjt -d 'Splice-junction error threshold (bp). (`-sjt`)' -r
complete -c tama -n "__fish_tama_using_subcommand collapse" -l lde -d 'Local density error threshold. (`-lde`)' -r
complete -c tama -n "__fish_tama_using_subcommand collapse" -l rm -d 'Run mode. (`-rm`)' -r -f -a "original\t'Hold everything in memory (default)'
low_mem\t'Lower-memory mode'"
complete -c tama -n "__fish_tama_using_subcommand collapse" -l vc -d 'Variation coverage threshold (reads). (`-vc`)' -r
complete -c tama -n "__fish_tama_using_subcommand collapse" -s b -l bam -d 'Treat input as BAM instead of SAM. (`-b`)'
complete -c tama -n "__fish_tama_using_subcommand collapse" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand collapse" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c tama -n "__fish_tama_using_subcommand merge" -s f -l filelist -d 'File list describing the annotations to merge. (`-f`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand merge" -s p -l prefix -d 'Output prefix. (`-p`)' -r
complete -c tama -n "__fish_tama_using_subcommand merge" -s e -l ends -d 'Collapse exon ends. (`-e`)' -r -f -a "common_ends\t'Use the most common start/end among the collapsed models'
longest_ends\t'Use the longest start/end among the collapsed models'"
complete -c tama -n "__fish_tama_using_subcommand merge" -s a -l five-prime -d '5\' threshold. (`-a`, merge default 20)' -r
complete -c tama -n "__fish_tama_using_subcommand merge" -s m -l exon-thresh -d 'Exon/splice-junction threshold. (`-m`)' -r
complete -c tama -n "__fish_tama_using_subcommand merge" -s z -l three-prime -d '3\' threshold. (`-z`, merge default 20)' -r
complete -c tama -n "__fish_tama_using_subcommand merge" -s d -l dup -d 'Duplicate merge behaviour. (`-d`)' -r -f -a "merge_dup\t'Merge duplicate models into one'
no_merge\t'Keep duplicate models separate'"
complete -c tama -n "__fish_tama_using_subcommand merge" -s s -l source-id -d 'Use gene/transcript IDs from this merge source. (`-s`)' -r
complete -c tama -n "__fish_tama_using_subcommand merge" -l cds -d 'Use CDS from this merge source. (`-cds`)' -r
complete -c tama -n "__fish_tama_using_subcommand merge" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand merge" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c tama -n "__fish_tama_using_subcommand format; and not __fish_seen_subcommand_from bed2gtf bed2gtf-orf fastq2fasta gtf2bed gff2bed id-filter help" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand format; and not __fish_seen_subcommand_from bed2gtf bed2gtf-orf fastq2fasta gtf2bed gff2bed id-filter help" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand format; and not __fish_seen_subcommand_from bed2gtf bed2gtf-orf fastq2fasta gtf2bed gff2bed id-filter help" -f -a "bed2gtf" -d 'Convert TAMA BED to Ensembl-style GTF (no CDS). (tama_convert_bed_gtf_ensembl_no_cds)'
complete -c tama -n "__fish_tama_using_subcommand format; and not __fish_seen_subcommand_from bed2gtf bed2gtf-orf fastq2fasta gtf2bed gff2bed id-filter help" -f -a "bed2gtf-orf" -d 'Convert TAMA BED to GTF with ORF/NMD CDS. (tama_convert_bed_gtf_ensembl_orf_nmd)'
complete -c tama -n "__fish_tama_using_subcommand format; and not __fish_seen_subcommand_from bed2gtf bed2gtf-orf fastq2fasta gtf2bed gff2bed id-filter help" -f -a "fastq2fasta" -d 'Convert Nanopore FASTQ to FASTA. (tama_convert_nanopore_fastq_fasta)'
complete -c tama -n "__fish_tama_using_subcommand format; and not __fish_seen_subcommand_from bed2gtf bed2gtf-orf fastq2fasta gtf2bed gff2bed id-filter help" -f -a "gtf2bed" -d 'Convert a GTF to TAMA BED12. (tama_format_gtf_to_bed12_*)'
complete -c tama -n "__fish_tama_using_subcommand format; and not __fish_seen_subcommand_from bed2gtf bed2gtf-orf fastq2fasta gtf2bed gff2bed id-filter help" -f -a "gff2bed" -d 'Convert a GFF to TAMA BED12. (tama_format_gff_to_bed12_*)'
complete -c tama -n "__fish_tama_using_subcommand format; and not __fish_seen_subcommand_from bed2gtf bed2gtf-orf fastq2fasta gtf2bed gff2bed id-filter help" -f -a "id-filter" -d 'Restructure/filter BED ID fields. (tama_format_id_filter)'
complete -c tama -n "__fish_tama_using_subcommand format; and not __fish_seen_subcommand_from bed2gtf bed2gtf-orf fastq2fasta gtf2bed gff2bed id-filter help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from bed2gtf" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from bed2gtf" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from bed2gtf-orf" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from bed2gtf-orf" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from fastq2fasta" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from fastq2fasta" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from gtf2bed" -l source -r -f -a "ensembl\t''
ncbi\t''
stringtie\t''"
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from gtf2bed" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from gtf2bed" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from gff2bed" -l source -r -f -a "cupcake\t''
liftoff\t''"
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from gff2bed" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from gff2bed" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from id-filter" -s b -l bed -d 'BED file. (`-b`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from id-filter" -s o -l output -d 'Output file. (`-o`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from id-filter" -s f -l filter -d 'Filter level. (`-f`)' -r -f -a "none\t'Keep every model, filling missing Ensembl IDs from the TAMA IDs'
only_match\t'Keep only models that carry an Ensembl ID'"
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from id-filter" -s s -l method -d 'Sub-field method. (`-s`)' -r -f -a "ensembl_merge\t'ID line produced by `tama merge` against an Ensembl annotation'
ensembl_orf\t'ID line produced by the ORF/NMD pipeline'
custom\t'Reorder arbitrary sub-fields with `--reshuffle`'"
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from id-filter" -s r -l reshuffle -d 'Custom reshuffle parameter, e.g. `3,4,1,2`. (`-r`)' -r
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from id-filter" -s d -l delim -d 'Sub-field delimiters. (`-d`)' -r
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from id-filter" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from id-filter" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from help" -f -a "bed2gtf" -d 'Convert TAMA BED to Ensembl-style GTF (no CDS). (tama_convert_bed_gtf_ensembl_no_cds)'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from help" -f -a "bed2gtf-orf" -d 'Convert TAMA BED to GTF with ORF/NMD CDS. (tama_convert_bed_gtf_ensembl_orf_nmd)'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from help" -f -a "fastq2fasta" -d 'Convert Nanopore FASTQ to FASTA. (tama_convert_nanopore_fastq_fasta)'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from help" -f -a "gtf2bed" -d 'Convert a GTF to TAMA BED12. (tama_format_gtf_to_bed12_*)'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from help" -f -a "gff2bed" -d 'Convert a GFF to TAMA BED12. (tama_format_gff_to_bed12_*)'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from help" -f -a "id-filter" -d 'Restructure/filter BED ID fields. (tama_format_id_filter)'
complete -c tama -n "__fish_tama_using_subcommand format; and __fish_seen_subcommand_from help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand filter; and not __fish_seen_subcommand_from primary-orf fragments polya single-read help" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand filter; and not __fish_seen_subcommand_from primary-orf fragments polya single-read help" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand filter; and not __fish_seen_subcommand_from primary-orf fragments polya single-read help" -f -a "primary-orf" -d 'Keep primary transcripts by ORF. (tama_filter_primary_transcripts_orf)'
complete -c tama -n "__fish_tama_using_subcommand filter; and not __fish_seen_subcommand_from primary-orf fragments polya single-read help" -f -a "fragments" -d 'Remove fragment models. (tama_remove_fragment_models)'
complete -c tama -n "__fish_tama_using_subcommand filter; and not __fish_seen_subcommand_from primary-orf fragments polya single-read help" -f -a "polya" -d 'Remove poly-A models by level. (tama_remove_polya_models_levels)'
complete -c tama -n "__fish_tama_using_subcommand filter; and not __fish_seen_subcommand_from primary-orf fragments polya single-read help" -f -a "single-read" -d 'Remove single-read models by level. (tama_remove_single_read_models_levels)'
complete -c tama -n "__fish_tama_using_subcommand filter; and not __fish_seen_subcommand_from primary-orf fragments polya single-read help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from primary-orf" -s b -l bed -d 'ORF/NMD BED file. (`-b`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from primary-orf" -s o -l output -d 'Output BED file. (`-o`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from primary-orf" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from primary-orf" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from fragments" -s f -l bed -d 'BED file. (`-f`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from fragments" -s o -l output -d 'Output prefix. (`-o`)' -r
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from fragments" -s m -l wobble -d 'Exon/splice-junction wobble threshold. (`-m`)' -r
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from fragments" -s e -l ends-wobble -d 'Transcript ends wobble threshold. (`-e`)' -r
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from fragments" -s s -l overlap-percent -d 'Single-exon overlap percent threshold. (`-s`)' -r
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from fragments" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from fragments" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from polya" -s b -l bed -d 'Annotation BED file. (`-b`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from polya" -s f -l filelist -d 'Filelist: `source_name<TAB>polya_file`. (`-f`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from polya" -s r -l read -d 'Read support levels file. (`-r`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from polya" -s o -l output -d 'Output prefix. (`-o`)' -r
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from polya" -s p -l percent -d 'Percent poly-A threshold. (`-p`)' -r
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from polya" -s l -l level -d 'Removal level. (`-l`)' -r -f -a "gene\t'Remove at the gene level'
transcript\t'Remove at the transcript level'"
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from polya" -s a -l support -d 'Which poly-A supported models to consider. (`-a`)' -r -f -a "all_polya\t'Consider every poly-A supported model'
singleton_polya\t'Consider only single-read poly-A supported models'"
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from polya" -s k -l multi -d 'Multi-exon handling. (`-k`)' -r -f -a "keep_multi\t'Keep multi-exon models regardless of support'
remove_multi\t'Subject multi-exon models to the same filter'"
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from polya" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from polya" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from single-read" -s b -l bed -d 'Annotation BED file. (`-b`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from single-read" -s r -l read -d 'Read support levels file. (`-r`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from single-read" -s o -l output -d 'Output prefix. (`-o`)' -r
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from single-read" -s l -l level -d 'Removal level. (`-l`)' -r -f -a "gene\t'Remove at the gene level'
transcript\t'Remove at the transcript level'"
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from single-read" -s k -l multi -d 'Multi-exon handling. (`-k`)' -r -f -a "keep_multi\t'Keep multi-exon models regardless of support'
remove_multi\t'Subject multi-exon models to the same filter'"
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from single-read" -s s -l source-support -d 'Minimum number of supporting sources. (`-s`)' -r
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from single-read" -s n -l read-support -d 'Minimum number of supporting reads. (`-n`)' -r
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from single-read" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from single-read" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from help" -f -a "primary-orf" -d 'Keep primary transcripts by ORF. (tama_filter_primary_transcripts_orf)'
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from help" -f -a "fragments" -d 'Remove fragment models. (tama_remove_fragment_models)'
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from help" -f -a "polya" -d 'Remove poly-A models by level. (tama_remove_polya_models_levels)'
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from help" -f -a "single-read" -d 'Remove single-read models by level. (tama_remove_single_read_models_levels)'
complete -c tama -n "__fish_tama_using_subcommand filter; and __fish_seen_subcommand_from help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand orf; and not __fish_seen_subcommand_from seek extract-cds add-cds blastp-parse help" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand orf; and not __fish_seen_subcommand_from seek extract-cds add-cds blastp-parse help" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand orf; and not __fish_seen_subcommand_from seek extract-cds add-cds blastp-parse help" -f -a "seek" -d 'Find ORFs in transcript sequences. (tama_orf_seeker)'
complete -c tama -n "__fish_tama_using_subcommand orf; and not __fish_seen_subcommand_from seek extract-cds add-cds blastp-parse help" -f -a "extract-cds" -d 'Extract CDS regions from a BED. (tama_bed_extract_cds)'
complete -c tama -n "__fish_tama_using_subcommand orf; and not __fish_seen_subcommand_from seek extract-cds add-cds blastp-parse help" -f -a "add-cds" -d 'Add CDS regions to a BED. (tama_cds_regions_bed_add)'
complete -c tama -n "__fish_tama_using_subcommand orf; and not __fish_seen_subcommand_from seek extract-cds add-cds blastp-parse help" -f -a "blastp-parse" -d 'Parse blastp output for ORF selection. (tama_orf_blastp_parser)'
complete -c tama -n "__fish_tama_using_subcommand orf; and not __fish_seen_subcommand_from seek extract-cds add-cds blastp-parse help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from seek" -s f -l fasta -d 'Transcript FASTA file. (`-f`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from seek" -s o -l output -d 'Output protein FASTA. (`-o`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from seek" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from seek" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from extract-cds" -s b -l bed -d 'ORF/NMD BED file. (`-b`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from extract-cds" -s s -l stop -d 'Whether the stop codon is part of the CDS. (`-s`)' -r -f -a "include_stop\t'Include the stop codon in the CDS'
no_stop_codon\t'Exclude the stop codon from the CDS'"
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from extract-cds" -s o -l output -d 'Output BED file. (`-o`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from extract-cds" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from extract-cds" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from add-cds" -s p -l parse -d 'Blastp parse file. (`-p`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from add-cds" -s a -l bed -d 'Annotation BED file. (`-a`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from add-cds" -s f -l fasta -d 'Transcript FASTA. (`-f`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from add-cds" -s o -l output -d 'Output BED file. (`-o`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from add-cds" -s s -l stop -d 'Whether the stop codon is part of the CDS. (`-s`)' -r -f -a "include_stop\t'Include the stop codon in the CDS'
no_stop_codon\t'Exclude the stop codon from the CDS'"
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from add-cds" -s d -l sj-dist -d 'Distance from last SJ to call NMD. (`-d`)' -r
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from add-cds" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from add-cds" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from blastp-parse" -s b -l blastp -d 'BLASTP pairwise output file. (`-b`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from blastp-parse" -s o -l output -d 'Output file. (`-o`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from blastp-parse" -s f -l format -d 'DB ID format. (`-f`)' -r -f -a "uniref\t'UniRef-style subject IDs'
ensembl\t'Ensembl-style subject IDs (`gene:` / `transcript:` fields)'"
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from blastp-parse" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from blastp-parse" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from help" -f -a "seek" -d 'Find ORFs in transcript sequences. (tama_orf_seeker)'
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from help" -f -a "extract-cds" -d 'Extract CDS regions from a BED. (tama_bed_extract_cds)'
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from help" -f -a "add-cds" -d 'Add CDS regions to a BED. (tama_cds_regions_bed_add)'
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from help" -f -a "blastp-parse" -d 'Parse blastp output for ORF selection. (tama_orf_blastp_parser)'
complete -c tama -n "__fish_tama_using_subcommand orf; and __fish_seen_subcommand_from help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand support; and not __fish_seen_subcommand_from collapse-cluster levels merge-collapse help" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand support; and not __fish_seen_subcommand_from collapse-cluster levels merge-collapse help" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand support; and not __fish_seen_subcommand_from collapse-cluster levels merge-collapse help" -f -a "collapse-cluster" -d 'Read support from a collapse trans_read.bed + cluster file. (tama_read_support_collapse_cluster)'
complete -c tama -n "__fish_tama_using_subcommand support; and not __fish_seen_subcommand_from collapse-cluster levels merge-collapse help" -f -a "levels" -d 'Read support levels file. (tama_read_support_levels)'
complete -c tama -n "__fish_tama_using_subcommand support; and not __fish_seen_subcommand_from collapse-cluster levels merge-collapse help" -f -a "merge-collapse" -d 'Read support after merging collapse outputs. (tama_read_support_merge_collapse)'
complete -c tama -n "__fish_tama_using_subcommand support; and not __fish_seen_subcommand_from collapse-cluster levels merge-collapse help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from collapse-cluster" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from collapse-cluster" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from levels" -s f -l filelist -d 'Filelist: `source_name<TAB>transread_file<TAB>file_type`. (`-f`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from levels" -s m -l merge -d 'Merge `_merge.txt`, or `no_merge`. (`-m`)' -r
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from levels" -s o -l output -d 'Output prefix (writes `<prefix>_read_support.txt`). (`-o`)' -r
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from levels" -l mt -d 'Merge file layout. (`-mt`)' -r -f -a "tama\t'A `tama merge` `_merge.txt` / `trans_read.bed` file'
cupcake\t'A cupcake collapse group file'
filter\t'A `tama filter` report'"
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from levels" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from levels" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from merge-collapse" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from merge-collapse" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from help" -f -a "collapse-cluster" -d 'Read support from a collapse trans_read.bed + cluster file. (tama_read_support_collapse_cluster)'
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from help" -f -a "levels" -d 'Read support levels file. (tama_read_support_levels)'
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from help" -f -a "merge-collapse" -d 'Read support after merging collapse outputs. (tama_read_support_merge_collapse)'
complete -c tama -n "__fish_tama_using_subcommand support; and __fish_seen_subcommand_from help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand stats; and not __fish_seen_subcommand_from degradation model-changes saturation help" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand stats; and not __fish_seen_subcommand_from degradation model-changes saturation help" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand stats; and not __fish_seen_subcommand_from degradation model-changes saturation help" -f -a "degradation" -d 'Degradation signature from capped vs no_cap collapse. (tama_degradation_signature)'
complete -c tama -n "__fish_tama_using_subcommand stats; and not __fish_seen_subcommand_from degradation model-changes saturation help" -f -a "model-changes" -d 'Find model changes between two sources. (tama_find_model_changes)'
complete -c tama -n "__fish_tama_using_subcommand stats; and not __fish_seen_subcommand_from degradation model-changes saturation help" -f -a "saturation" -d 'Sampling saturation curve from a read-support levels file. (tama_sampling_saturation_curve)'
complete -c tama -n "__fish_tama_using_subcommand stats; and not __fish_seen_subcommand_from degradation model-changes saturation help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from degradation" -s c -l capped -d 'Capped collapse `trans_read.bed`. (`-c`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from degradation" -l nc -d 'No-cap collapse `trans_read.bed`. (`-nc`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from degradation" -s o -l output -d 'Output file. (`-o`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from degradation" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from degradation" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from model-changes" -s b -l bed -d 'Annotation BED file. (`-b`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from model-changes" -s r -l read -d 'Read support levels file. (`-r`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from model-changes" -s o -l output -d 'Output prefix. (`-o`)' -r
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from model-changes" -l ref -d 'Reference source name. (`-ref`)' -r
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from model-changes" -l alt -d 'Alternative source name. (`-alt`)' -r
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from model-changes" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from model-changes" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from saturation" -s r -l report -d 'Read support levels file. (`-r`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from saturation" -s b -l bin -d 'Read bin size. (`-b`)' -r
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from saturation" -s o -l output -d 'Output file. (`-o`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from saturation" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from saturation" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from help" -f -a "degradation" -d 'Degradation signature from capped vs no_cap collapse. (tama_degradation_signature)'
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from help" -f -a "model-changes" -d 'Find model changes between two sources. (tama_find_model_changes)'
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from help" -f -a "saturation" -d 'Sampling saturation curve from a read-support levels file. (tama_sampling_saturation_curve)'
complete -c tama -n "__fish_tama_using_subcommand stats; and __fish_seen_subcommand_from help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand split; and not __fish_seen_subcommand_from fasta sam help" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand split; and not __fish_seen_subcommand_from fasta sam help" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand split; and not __fish_seen_subcommand_from fasta sam help" -f -a "fasta" -d 'Split a FASTA into chunks. (tama_fasta_splitter)'
complete -c tama -n "__fish_tama_using_subcommand split; and not __fish_seen_subcommand_from fasta sam help" -f -a "sam" -d 'Split a mapped SAM by chromosome. (tama_mapped_sam_splitter)'
complete -c tama -n "__fish_tama_using_subcommand split; and not __fish_seen_subcommand_from fasta sam help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand split; and __fish_seen_subcommand_from fasta" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand split; and __fish_seen_subcommand_from fasta" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand split; and __fish_seen_subcommand_from sam" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand split; and __fish_seen_subcommand_from sam" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand split; and __fish_seen_subcommand_from help" -f -a "fasta" -d 'Split a FASTA into chunks. (tama_fasta_splitter)'
complete -c tama -n "__fish_tama_using_subcommand split; and __fish_seen_subcommand_from help" -f -a "sam" -d 'Split a mapped SAM by chromosome. (tama_mapped_sam_splitter)'
complete -c tama -n "__fish_tama_using_subcommand split; and __fish_seen_subcommand_from help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand variants; and not __fish_seen_subcommand_from call help" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand variants; and not __fish_seen_subcommand_from call help" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand variants; and not __fish_seen_subcommand_from call help" -f -a "call" -d 'Call variants from a sorted SAM against the genome. (tama_variant_caller)'
complete -c tama -n "__fish_tama_using_subcommand variants; and not __fish_seen_subcommand_from call help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from call" -s s -l sam -d 'Sorted SAM file. (`-s`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from call" -s f -l fasta -d 'Genome FASTA file. (`-f`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from call" -s p -l prefix -d 'Output prefix. (`-p`)' -r
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from call" -s x -l cap-flag -d 'Capped flag. (`-x`)' -r -f -a "capped\t'5\' cap-trimmed data; 5\' ends are trusted'
no_cap\t'Non-cap-trimmed data; 5\' degradation is expected'"
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from call" -s c -l coverage -d 'Minimum coverage percent. (`-c`)' -r
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from call" -s i -l identity -d 'Minimum identity percent. (`-i`)' -r
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from call" -l icm -d 'Identity method. (`-icm`)' -r -f -a "ident_cov\t'Include hard/soft clipping in the denominator (default)'
ident_map\t'Exclude hard/soft clipping'"
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from call" -l sjt -d 'Splice-junction error threshold (bp). (`-sjt`)' -r
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from call" -l vc -d 'Variation coverage threshold (reads). (`-vc`)' -r
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from call" -s b -l bam -d 'Treat input as BAM instead of SAM. (`-b`)'
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from call" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from call" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from help" -f -a "call" -d 'Call variants from a sorted SAM against the genome. (tama_variant_caller)'
complete -c tama -n "__fish_tama_using_subcommand variants; and __fish_seen_subcommand_from help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand cleanup; and not __fish_seen_subcommand_from polya help" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand cleanup; and not __fish_seen_subcommand_from polya help" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand cleanup; and not __fish_seen_subcommand_from polya help" -f -a "polya" -d 'Trim poly-A tails from FLNC reads. (tama_flnc_polya_cleanup)'
complete -c tama -n "__fish_tama_using_subcommand cleanup; and not __fish_seen_subcommand_from polya help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand cleanup; and __fish_seen_subcommand_from polya" -s f -l fasta -d 'FLNC FASTA file. (`-f`)' -r -F
complete -c tama -n "__fish_tama_using_subcommand cleanup; and __fish_seen_subcommand_from polya" -s p -l prefix -d 'Output prefix. (`-p`)' -r
complete -c tama -n "__fish_tama_using_subcommand cleanup; and __fish_seen_subcommand_from polya" -s m -l min-length -d 'Minimum read length to keep. (`-m`)' -r
complete -c tama -n "__fish_tama_using_subcommand cleanup; and __fish_seen_subcommand_from polya" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand cleanup; and __fish_seen_subcommand_from polya" -s h -l help -d 'Print help'
complete -c tama -n "__fish_tama_using_subcommand cleanup; and __fish_seen_subcommand_from help" -f -a "polya" -d 'Trim poly-A tails from FLNC reads. (tama_flnc_polya_cleanup)'
complete -c tama -n "__fish_tama_using_subcommand cleanup; and __fish_seen_subcommand_from help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand completions" -s v -l verbose -d 'Print progress messages to stderr while running (so you can see the run is still going). Off by default; only end-of-run summaries are shown. Can also be controlled with the RUST_LOG environment variable'
complete -c tama -n "__fish_tama_using_subcommand completions" -s h -l help -d 'Print help (see more with \'--help\')'
complete -c tama -n "__fish_tama_using_subcommand help; and not __fish_seen_subcommand_from collapse merge format filter orf support stats split variants cleanup completions help" -f -a "collapse" -d 'Collapse mapped reads into transcript models (tama_collapse)'
complete -c tama -n "__fish_tama_using_subcommand help; and not __fish_seen_subcommand_from collapse merge format filter orf support stats split variants cleanup completions help" -f -a "merge" -d 'Merge transcript annotations across sources (tama_merge)'
complete -c tama -n "__fish_tama_using_subcommand help; and not __fish_seen_subcommand_from collapse merge format filter orf support stats split variants cleanup completions help" -f -a "format" -d 'Format converters (BED <-> GTF/GFF, FASTQ -> FASTA, ...)'
complete -c tama -n "__fish_tama_using_subcommand help; and not __fish_seen_subcommand_from collapse merge format filter orf support stats split variants cleanup completions help" -f -a "filter" -d 'Filter transcript models'
complete -c tama -n "__fish_tama_using_subcommand help; and not __fish_seen_subcommand_from collapse merge format filter orf support stats split variants cleanup completions help" -f -a "orf" -d 'ORF / NMD prediction tools'
complete -c tama -n "__fish_tama_using_subcommand help; and not __fish_seen_subcommand_from collapse merge format filter orf support stats split variants cleanup completions help" -f -a "support" -d 'Read-support tracking tools'
complete -c tama -n "__fish_tama_using_subcommand help; and not __fish_seen_subcommand_from collapse merge format filter orf support stats split variants cleanup completions help" -f -a "stats" -d 'File statistics tools'
complete -c tama -n "__fish_tama_using_subcommand help; and not __fish_seen_subcommand_from collapse merge format filter orf support stats split variants cleanup completions help" -f -a "split" -d 'File splitting tools'
complete -c tama -n "__fish_tama_using_subcommand help; and not __fish_seen_subcommand_from collapse merge format filter orf support stats split variants cleanup completions help" -f -a "variants" -d 'Variant calling'
complete -c tama -n "__fish_tama_using_subcommand help; and not __fish_seen_subcommand_from collapse merge format filter orf support stats split variants cleanup completions help" -f -a "cleanup" -d 'Sequence cleanup tools'
complete -c tama -n "__fish_tama_using_subcommand help; and not __fish_seen_subcommand_from collapse merge format filter orf support stats split variants cleanup completions help" -f -a "completions" -d 'Generate a shell completion script on stdout'
complete -c tama -n "__fish_tama_using_subcommand help; and not __fish_seen_subcommand_from collapse merge format filter orf support stats split variants cleanup completions help" -f -a "help" -d 'Print this message or the help of the given subcommand(s)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from format" -f -a "bed2gtf" -d 'Convert TAMA BED to Ensembl-style GTF (no CDS). (tama_convert_bed_gtf_ensembl_no_cds)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from format" -f -a "bed2gtf-orf" -d 'Convert TAMA BED to GTF with ORF/NMD CDS. (tama_convert_bed_gtf_ensembl_orf_nmd)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from format" -f -a "fastq2fasta" -d 'Convert Nanopore FASTQ to FASTA. (tama_convert_nanopore_fastq_fasta)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from format" -f -a "gtf2bed" -d 'Convert a GTF to TAMA BED12. (tama_format_gtf_to_bed12_*)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from format" -f -a "gff2bed" -d 'Convert a GFF to TAMA BED12. (tama_format_gff_to_bed12_*)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from format" -f -a "id-filter" -d 'Restructure/filter BED ID fields. (tama_format_id_filter)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from filter" -f -a "primary-orf" -d 'Keep primary transcripts by ORF. (tama_filter_primary_transcripts_orf)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from filter" -f -a "fragments" -d 'Remove fragment models. (tama_remove_fragment_models)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from filter" -f -a "polya" -d 'Remove poly-A models by level. (tama_remove_polya_models_levels)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from filter" -f -a "single-read" -d 'Remove single-read models by level. (tama_remove_single_read_models_levels)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from orf" -f -a "seek" -d 'Find ORFs in transcript sequences. (tama_orf_seeker)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from orf" -f -a "extract-cds" -d 'Extract CDS regions from a BED. (tama_bed_extract_cds)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from orf" -f -a "add-cds" -d 'Add CDS regions to a BED. (tama_cds_regions_bed_add)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from orf" -f -a "blastp-parse" -d 'Parse blastp output for ORF selection. (tama_orf_blastp_parser)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from support" -f -a "collapse-cluster" -d 'Read support from a collapse trans_read.bed + cluster file. (tama_read_support_collapse_cluster)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from support" -f -a "levels" -d 'Read support levels file. (tama_read_support_levels)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from support" -f -a "merge-collapse" -d 'Read support after merging collapse outputs. (tama_read_support_merge_collapse)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from stats" -f -a "degradation" -d 'Degradation signature from capped vs no_cap collapse. (tama_degradation_signature)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from stats" -f -a "model-changes" -d 'Find model changes between two sources. (tama_find_model_changes)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from stats" -f -a "saturation" -d 'Sampling saturation curve from a read-support levels file. (tama_sampling_saturation_curve)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from split" -f -a "fasta" -d 'Split a FASTA into chunks. (tama_fasta_splitter)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from split" -f -a "sam" -d 'Split a mapped SAM by chromosome. (tama_mapped_sam_splitter)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from variants" -f -a "call" -d 'Call variants from a sorted SAM against the genome. (tama_variant_caller)'
complete -c tama -n "__fish_tama_using_subcommand help; and __fish_seen_subcommand_from cleanup" -f -a "polya" -d 'Trim poly-A tails from FLNC reads. (tama_flnc_polya_cleanup)'
complete -c tama -n "__fish_tama_using_subcommand completions" -f -a "bash elvish fish powershell zsh"
