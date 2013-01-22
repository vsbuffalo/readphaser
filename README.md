# readphaser -- separate reads based on mapping results and hapcut data

`phasereads.py` takes phased variants from HapCut and a corresponding
BAM file of alignments to produce a FASTA file of phased reads (with
phasing block and haplotype in the header).

## Implementation 

## Limitations

**Indels** are not handled: `phasereads.py` assumes a phased variant
  has the same width in reads as it does on the reference. Multibase
  variants (i.e. from cluttered SNPs) will be handled as these are the
  same width in reads as in the reference. However, reads that
  partially overlap a mulitbase variant are not used in
  `phasereads.py` as they are ambigious. 

Indels will likely be implicitly handled in the event they are part of
reads that contain a phased variant. After phasing, reads containing
this phased variant may be used to regenerate a phased sequence
including the variant.

Indels are also a complexity around SNPs and multibase variants. In
some cases, a biallelic multibase can be close to an indel (due to
alignment issues, it may be higher scoring to open a gap than allow a
single mismatch).

## Statistics Output

To ensure reads are being grouped properly, a file of the 

