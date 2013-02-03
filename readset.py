"""
readset.py contains one class: ReadSet, which is a small abstraction
of sets to reads, helping with read accounting and ensuring reads
encounted twice during read separation are not output twice.

"""

from Bio.Seq import Seq
from collections import namedtuple

Read = namedtuple("Read", ["seq", "qual"])

def revcomp(seq):
    return str(Seq(seq).reverse_complement())

class ReadSet(object):
    """
    ReadSets are a lightweight way of storing paired-end reads in a
    way that ensures the same read (identified by read name) is not
    stored more than once (thus, a set). This encapsulates a dict with
    read names as the key and a list of two as the value. The list is
    to store pairs, in order. Each pair is composed of a tuple of
    (seq, qual).
    """
    def __init__(self, **kwargs):
        self.readset = dict()
        self._keyvals = kwargs

    def add_keyval(self, **kwargs):
        self._keyvals = dict(self._keyvals.items() + kwargs.items())

    @property
    def keyvals(self):
        return " ".join("%s:%s" % (k, v) for k, v in self._keyvals.items())
    
    def add_readpair(self, readpair):
        """
        Add a readpair (tuple of AlignedRead objects from pysam).
        """
        qnames = list(set(read.qname for read in readpair if read is not None))
        assert(len(qnames) == 1)
        qname = qnames[0]
        if self.readset.get(qname, None) is None:
            self.readset[qname] = [None, None]
        for read in readpair:
            if read is None:
                continue
            seq = read.seq
            qual = read.qual
            if read.is_reverse:
                seq = revcomp(seq)
                qual = qual[::-1]
            which_read = 0 if read.is_read1 else 1
            if None not in self.readset[read.qname]:
                # check that if there's a sequence in there, it's the same
                # sequence
                assert(self.readset[read.qname][which_read].seq == seq)
            else:
                self.readset[read.qname][which_read] = Read(seq, qual)

    def __len__(self):
        return len(self.readset)

    def __str__(self):
        """
        FASTQ stringified version of reads.
        """
        out = list()
        for qname, readpair in self.readset.iteritems():
            for i, read in enumerate(readpair):
                if read is None:
                    continue
                which = i + 1
                header = "%s-%s %s" % (qname, which, self.keyvals)
                out.append("@%s\n%s\n+\n%s\n" % (header, read.seq, read.qual))
        return out
                
    def write(self, file_handle):
        """
        Write FASTQ reads to file.
        """
        for qname, readpair in self.readset.iteritems():
            for i, read in enumerate(readpair):
                if read is None:
                    continue
                which = i + 1
                header = "%s-%s %s" % (qname, which, self.keyvals)
                file_handle.write("@%s\n%s\n+\n%s\n" % (header, read.seq, read.qual))

    def __iter__(self):
        """
        Iterate through a readset, yielding each read. This is
        irrespective of which member of a pair it is.
        """
        for readpair in self.readset.values():
            # read_group is both pairs
            for read in readpair:
                if read is not None:
                    yield read
