"""
pr.py -- phasereads.py

hapcut.HapCut() takes a hapcut file and returns a dictionary of
refnames and a list of hapcut.Blocks. Each hapcut.Block contains
block info and a list of block entries (phased variants).
"""
TEST = True
import pdb
import sys
import argparse
from collections import defaultdict, OrderedDict, Counter
from operator import itemgetter
import pysam
from Bio.Seq import Seq
from hapcut import HapCut
import fermi as fm

def revcomp(seq):
    return str(Seq(seq).reverse_complement())

class ReadSet(object):
    def __init__(self, refname, block, haplotype):
        self.refname = refname
        self.block = block
        self.haplotype = haplotype
        self.readset = dict()
    def add_read(self, read):
        """
        Add a read, modifying header to include the contig and phasing
        info. Also, look at mapping orientation and reverse if necessary.
        """
        if read.qname not in self.readset:
            self.readset[read.qname] = [None, None]
        which_read = 0 if read.is_read1 else 1
        seq = read.seq
        qual = read.qual
        if read.is_reverse:
            seq = revcomp(seq)
            qual = qual[::-1]
        if TEST and None not in self.readset[read.qname]:
            assert(self.readset[read.qname][which_read][0] == seq)
        self.readset[read.qname][which_read] = (seq, qual)

    @property
    def name(self):
        fields = (self.refname, self.block, self.haplotype)
        return "CT:%s BL:%s PH:%s" % fields
    
    def write_reads(self, file_handle):
        for qname, seqs in self.readset.iteritems():
            for i, seq in enumerate(seqs):
                if seq is None:
                    continue
                which = i + 1
                fields = map(str, (qname, which, self.refname, self.block, self.haplotype))
                header = "%s-%s CT:%s BL:%s PH:%s" % tuple(fields)
                file_handle.write("@%s\n%s\n+\n%s\n" % (header, seq[0], seq[1]))
    def __iter__(self):
        for read_group in self.readset.values():
            # read_group is both pairs
            for read in read_group:
                if read is not None:
                    yield read

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
    for read in bamfile.fetch(reference=refname):
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

def print_block_stats(refname, block_id, allele_counts, read_stats):
    """
    For each block, print allele counts and the read statistics.
    """
    tmp = (refname, block_id, read_stats['phased'], read_stats['indel_ignore'], read_stats['unused'])
    sys.stdout.write("# contig='%s' block_id=%d num_phased=%d indel_ignore=%d unused=%d\n" % tmp)
    sorted_counts = sorted(allele_counts.items(), key=itemgetter(0))
    for pos, counts in sorted_counts:
        joined = ";".join(["%s:%s" % (a, c) for a, c in counts.items()])
        sys.stdout.write("%s\t%d\t%d\t%s\n" % (refname, block_id, pos, joined))
    sys.stdout.flush()
    
def group_reads_by_block(block, block_id, bamfile, mapq, exclude_duplicates, getmate, callback=None):
    """
    Take a HapCut hapcut.Block (a single phased block) and a BAM file
    and group the BAM file's reads.

    Intervals are 0-indexed. HapCut is 1-indexed (thus,
    entry.position-1). 
    """
    refname = block.entries[0].chromosome
    sys.stderr.write("[phase_reads] phasing '%s', block_id %d\n" % (refname, block_id))
    sys.stderr.flush()
    # make a dictionary of all variants in a phased block
    haplotypes = dict()
    for entry in block.entries:
        ref_key = (entry.ref_allele, entry.haplotype_1)
        var_key = (entry.var_allele, entry.haplotype_2)
        allele_len = len(entry.ref_allele)
        haplotypes[(entry.position-1, entry.position-1 + allele_len, allele_len)] = dict([ref_key, var_key])

    readsets = (ReadSet(refname, block_id, 0), ReadSet(refname, block_id, 1))
    unused_phased = ReadSet(refname, "NA", "NA")
    stats = Counter()
    allele_counts = defaultdict(Counter)
                
    for read in bamfile.fetch(refname):
        if read.is_unmapped or read.mapq < mapq or (exclude_duplicates and read.is_duplicate):
            continue
        for key, alleles in haplotypes.items():
            interval = key[0:2]
            allele_len = key[2]
            if read.overlap(interval[0], interval[1]) == allele_len:
                # full overlap of variant - necessary for multibase variants

                # check for indel, which will break our variant
                # retrieval.
                if has_indel(read):
                    # TODO: we are losing reads here... possibly many
                    unused_phased.add_read(read)
                    stats['indel_ignore'] += 1
                else:
                    # no indel; can safely grab variant with position offet
                    read_var = read.query[interval[0]-read.pos:interval[1]-read.pos]
                    htype = alleles.get(read_var, None)
                    if htype is not None:
                        stats['phased'] += 1
                        readsets[htype].add_read(read)
                        allele_counts[interval[0]][read_var] += 1
                        read_mate = getmate(read)
                        if read_mate is not None:
                            readsets[htype].add_read(read_mate)
                            stats['phased_mates'] += 1
                    else:
                        # read does not have a phased variant
                        unused_phased.add_read(read)
                        stats['unused'] += 1
    print_block_stats(refname, block_id, allele_counts, stats)
    # TODO how to handle mates for phased and unphased
    if callback is not None:
        callback(readsets, unused_phased)

def phase_reads(bam_filename, hapcut_file, unphased_file, mapq, exclude_duplicates,
                callback, region=None):

    sys.stderr.write("[phase_reads] opening alignment BAM file...\t")
    bamfile = pysam.Samfile(bam_filename, 'rb')
    sys.stderr.write("done.\n")

    hapcut_dict = HapCut(hapcut_file).to_dict()

    if region is not None:
        hapcut_dict = dict([(region, hapcut_dict[region])])

    for refname, phased_blocks in hapcut_dict.iteritems():
        getmate = getmate_factory(refname, bamfile)
        for block_id, block in enumerate(phased_blocks):
            group_reads_by_block(block, block_id, bamfile, mapq, exclude_duplicates,
                                 getmate, callback)
    
    # handle unphased contigs
    if unphased_file is not None:
        unphased_contigs = set(bamfile.references) - set(hapcut_dict.keys())
        for read in bamfile:
            if (read not in unphased_contigs or read.is_unmapped or
                read.mapq < mapq or (exclude_duplicates and read.is_duplicate)):
                continue
            which_read = 1 if read.is_read1 else 2
            seq = read.seq
            qual = read.qual
            if read.is_reverse:
                seq = revcomp(seq)
                qual = qual[::-1]
            fields = map(str, (read.qname, which_read, bamfile.getrname[read.tid]))
            # NP: not phased
            header = "%s-%s CT:NP BL:NP PH:NP" % tuple(fields)
            unphased_file.write("@%s\n%s\n+\n%s\n" % (header, seq[0], seq[1]))

def assemble_main(args):
    """
    main() function for assembly with fermi. This also defines a
    closure callback function over some of the arguments.
    """
    
    def assembly_callback(phased_readset, unused_readset):
        if args.unused_phased is not None:
            unused_readset.write_reads(args.unused_phased)
        for readset in phased_readset:
            fermi = fm.Fermi()
            for seq, qual in readset:
                fermi.addseq(seq, qual)
            fermi.correct()
            tigs = fermi.assemble(do_clean=True)
            root_name = readset.name
            tigs = fermi.fastq_to_list(tigs, root_name)
            for tig in tigs:
                args.contigs.write(">%s\n%s\n" % (tig.header, tig.seq))

    phase_reads(args.bam, args.hapcut, args.unphased, args.mapq,
                args.exclude_duplicates, assembly_callback,
                region=args.region)

def output_main(args):
    """
    main() function for outputting reads to a file. This also defines
    a closure callback function over some of the arguments.
    """
    
    def writer_callback(phased_readset, unused_readset):
        if args.phased is not None:
            phased_readset[0].write_reads(args.phased)
            phased_readset[1].write_reads(args.phased)
        if args.unused_phased is not None:
            unused_readset.write_reads(args.unused_phased)

    phase_reads(args.bam, args.hapcut, args.unphased, args.mapq,
                args.exclude_duplicates, writer_callback,
                region=args.region)

if __name__ == "__main__":
    msg = "divide reads into groups, based on HapCut phasing results"
    parser = argparse.ArgumentParser(description=msg)
    subparsers = parser.add_subparsers()
    parser_assemble = subparsers.add_parser('assemble', help='assemble reads with fermi')
    parser_assemble.add_argument("-m", "--mapq",
                                 help="mapping quality threshold (exclude if below)", 
                                 type=int, required=False, default=0)
    parser_assemble.add_argument("-d", "--exclude-duplicates", help="exclude duplicate reads", 
                                 action="store_true", default=True)
    parser_assemble.add_argument("-c", "--contigs",
                                 help="FASTA filename for assembled phasedc contigs",
                                 type=argparse.FileType('w'), default=None)
    parser_assemble.add_argument("-o", "--unused-phased",
                                 help="FASTA filename for reads from phased "
                                 "contigs unused during phasing",
                                 type=argparse.FileType('w'), default=None)
    parser_assemble.add_argument("-u", "--unphased",
                                 help="FASTA filename for reads from unphased contigs",
                                 type=argparse.FileType('w'), default=None)

    parser_output = subparsers.add_parser('output', help='output phased reads to file')
    parser_output.add_argument("-u", "--unphased",
                               help="FASTA filename for reads from unphased contigs",
                               type=argparse.FileType('w'), default=None, required=False)
    parser_output.add_argument("-p", "--phased",
                               help="FASTA filename for reads from phased contigs",
                               type=argparse.FileType('w'), default=None)
    parser_output.add_argument("-o", "--unused-phased",
                               help="FASTA filename for reads from phased "
                               "contigs unused during phasing",
                               type=argparse.FileType('w'), default=None)
    parser_output.add_argument("-m", "--mapq",
                               help="mapping quality threshold (exclude if below)", 
                               type=int, required=False, default=0)
    parser_output.add_argument("-d", "--exclude-duplicates", help="exclude duplicate reads", 
                               action="store_true", default=True)
    parser.add_argument("hapcut", help="hapcut file", default=None,
                        type=argparse.FileType('r'))
    parser.add_argument("bam", help="BAM file of aligned reads", default=None,
                        type=str)
    parser.add_argument("region", help="optional region", default=None,
                        type=str, nargs="?")
    parser_output.set_defaults(func=output_main)
    parser_assemble.set_defaults(func=assemble_main)
    
    args = parser.parse_args()
    args.func(args)
