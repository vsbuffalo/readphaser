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

    def __getitem__(self, key):
        return self._keyvals[key]

    @property
    def keyvals(self):
        return " ".join("%s:%s" % (k, v) for k, v in self._keyvals.items())

    def join(self, other_readset):
        """
        Join two readsets.
        """
        for qname, readpair in other_readset.readset.items():
            value = self.readset.get(qname, None)
            if value is None:
                self.readset[qname] = readpair
            elif None in value:
                # merge the two
                for which_read in (0, 1):
                    if self.readset[qname][which_read] is not None:
                        if other_readset[qname][which_read] is not None:
                            assert(other_readset[qname][which_read].seq == self.readset[qname][which_read].seq)
                            assert(other_readset[qname][which_read].qual == self.readset[qname][which_read].qual)
                    else:
                        if other_readset[qname][which_read] is not None:
                            self.readset[qname][which_read] = other_readset[qname][which_read]
            else:
                pass # no merging to be done
    
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
            seq = read.query
            qual = read.qqual
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
                header = "%s/%i %s" % (qname, i + 1, self.keyvals)
                entry = "@%s\n%s\n+\n%s" % (header, read.seq, read.qual)
                out.append(entry)
        return "\n".join(out) + "\n"
                
    def write(self, file_handle):
        """
        Write FASTQ reads to file.
        """
        for qname, readpair in self.readset.iteritems():
            for i, read in enumerate(readpair):
                if read is None:
                    continue
                header = "%s/%i %s" % (qname, i + 1, self.keyvals)
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
