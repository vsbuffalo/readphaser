# readphaser -- separate reads based on mapping results and hapcut data

`pr.py` takes phased variants from HapCut and a corresponding
BAM file of alignments to produce:

 - FASTA file of phased reads (with phasing block and haplotype in the
    header).

 - FASTA file of reads unused in phasing (due to sequencing errors at
    the same location of a variant, for example).

 - FASTA file of reads from unphased contigs.

Optionally, `pr.py` can call the Fermi assembler to assembly reads
directly.

## Requirements

 - pysam

## Limitations

 - Indels are not handled: `pr.py` currently ignores reads with
  indels, because they displace other variants.


## Todo

 - Refactor to use two distinct subcommands: one for outputting reads,
   one for assembling reads via fermi.

 - Add asynchronous callback functions for assembly that could be used
   for running assemblies during phasing.

 - Multiprocessing.