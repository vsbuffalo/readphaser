"""
pr.py -- phasereads.py

hapcut.HapCut() takes a hapcut file and returns a dictionary of
refnames and a list of hapcut.Blocks. Each hapcut.Block contains
block info and a list of block entries (phased variants).
"""

import pdb
import sys
import argparse
from collections import defaultdict, OrderedDict, Counter
import pysam
from Bio.Seq import Seq
from hapcut import HapCut
import fermi as fm

def revcomp(seq):
    return str(Seq(seq).reverse_complement())

def has_indel(read):
    return any([op in (1, 2) for op, _ in read.cigar])

def getmate_factory(refname, bamfile):
    """
    Hash all reads from a refname and return a function that grabs the
    mate given a read name and which read is currently being handled.

    Note that we explicitly do not handle KeyErrors here; this is
    intentional; if a KeyError occurs it means the mate is not mapped,
    or is mapped to a different contig (improper pair), which should
    be handled elsewhere.
    """
    reads = defaultdict(dict)
    for read in bamfile.fetch(region=refname):
        which_read = 0 if read.is_read1 else 1
        reads[read.qname][which_read] = read

    def getter(read):
        which_read = 0 if read.is_read1 else 1
        assert(not read.is_unmapped)
        usable_mate = not read.mate_is_unmapped and read.rnext == read.tid
        if not usable_mate:
            return None
        which = (1, 0)
        return reads[read.qname][which[which_read]]

    return getter


def group_reads_by_block(block, bamfile, getmate):
    """
    Take a HapCut hapcut.Block (a single phased block) and a BAM file
    and group the BAM file's reads.

    Intervals are 0-indexed. HapCut is 1-indexed (thus,
    entry.position-1). 
    """
    # make a dictionary of all variants in a phased block
    haplotypes = dict()
    for entry in block.entries:
        ref_key = (entry.ref_allele, entry.haplotype_1)
        var_key = (entry.var_allele, entry.haplotype_2)
        allele_len = len(entry.ref_allele)
        haplotypes[(entry.position-1, entry.position-1 + allele_len, allele_len)] = dict([ref_key, var_key])

    readsets = {0:list(), 1:list()} # TODO make set
    unused_phased = list()
    stats = Counter()
    for read in bamfile.fetch(entry.chromosome):
        if read.is_unmapped:
            continue
        for key, alleles in haplotypes.items():
            interval = key[0:2]
            allele_len = key[2]
            if read.overlap(interval[0], interval[1]) == allele_len:
                # full overlap of variant - necessary for multibase variants

                # check for indel, which will break our variant
                # retrieval.
                if has_indel(read):
                    stats['indel_ignore'] += 1
                else:
                    read_var = read.query[interval[0]-read.pos:interval[1]-read.pos]
                    htype = alleles.get(read_var, None)
                    if htype is not None:
                        stats['phased'] += 1
                        readsets[htype].append(read)
                        read_mate = getmate(read)
                        if read_mate is not None:
                            readsets[htype].append(read_mate)
                            stats['phased_mates'] += 1
                    else:
                        # read does not have a phased variant
                        unused_phased.append(read)
                        stats['unused'] += 1
                        pdb.set_trace()

    fermi = fm.Fermi()
    for read in readsets[0]:
        fermi.addseq(read.query)
    fermi.correct()
    tigs = fermi.assemble(do_clean=True)
    pdb.set_trace()

    # TODO how to handle mates for phased and unphased
    

def phase_reads(bam_filename, hapcut_filename, region=None):
    sys.stderr.write("[phase_reads] opening alignment BAM file...\t")
    bamfile = pysam.Samfile(bam_filename, 'rb')
    sys.stderr.write("done.\n")

    hapcut_dict = HapCut(open(hapcut_filename, 'r')).to_dict()

    if region is not None:
        hapcut_dict = dict([(region, hapcut_dict[region])])

    for refname, phased_blocks in hapcut_dict.iteritems():
        getmate = getmate_factory(refname, bamfile)
        for i, block in enumerate(phased_blocks):
            group_reads_by_block(block, bamfile, getmate)
    
    pdb.set_trace()
    


if __name__ == "__main__":
    msg = "divide reads into groups, based on HapCut phasing results"
    parser = argparse.ArgumentParser(description=msg)
    parser.add_argument("-u", "--unphased", help="FASTA filename for reads from unphased contigs",
                        type=str, default=None, required=False)
    parser.add_argument("-p", "--phased", help="FASTA filename for reads from phased contigs",
                        type=str, required=True)
    parser.add_argument("-o", "--unused-phased",
                        help="FASTA filename for reads from phased contigs unused during phasing",
                        type=str, default=None, required=False)
    parser.add_argument("-m", "--mapq", help="mapping quality threshold (exclude if below)", 
                        type=int, required=False, default=None)
    parser.add_argument("-d", "--exclude-duplicates", help="exclude duplicate reads", 
                        action="store_true", default=True)
    parser.add_argument("hapcut", help="hapcut file", default=None,
                        type=str)
    parser.add_argument("bam", help="BAM file of aligned reads", default=None,
                        type=str)
    parser.add_argument("region", help="optional region", default=None,
                        type=str, nargs="?")

    args = parser.parse_args()
    phase_reads(args.bam, args.hapcut, region=args.region)
