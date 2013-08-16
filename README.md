# readphaser -- separate reads based on mapping results and HapCUT data

**Caveat emptor**: `readphaser` was tested for use in
[Krasileva et al, 2013](http://genomebiology.com/2013/14/6/R66/abstract)
and benchmarked using data from this project. In our tests, it phases
many contigs where there is a sufficient density of variants and
coverage. However, it has not been widely tested in other species and
your mileage may vary. As with all bioinformatics program, check that
your results are consistent with your own validation procedures. Also,
I am quite busy nowadays so (as with all free software) there is no
guarantee of support, but I will try my best.

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
