# readphaser -- separate reads based on mapping results and HapCUT data

`pr.py` takes phased variants from
[HapCut](https://sites.google.com/site/vibansal/software/hapcut) and a
corresponding BAM file of alignments to produce:

 - FASTA file of phased reads (with phasing block and haplotype in the
    header).

 - FASTA file of reads unused in phasing (due to sequencing errors at
    the same location of a variant, for example).

 - FASTA file of reads from unphased contigs.

## Requirements

 - pysam

## Limitations

 - Indels are not handled: `pr.py` currently ignores reads with
  indels, because (1) they displace other variants and (2) variants
  around indels are often incorrectly called.