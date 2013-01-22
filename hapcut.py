"""
hapcut.py - iterate through HapCut output

Descriptions of HapCut output from readme:

## Block Header
 - offset: first_variant_block
 - len: length_of_block
 - phased: phased_variants_block
 - SPAN: lengthspanned
 - MECscore score
 - fragments #fragments

## Block Entries 
 - variant_id
 - haplotype_1
 - haplotype_2
 - chromosome
 - position
 - refallele
 - variantallele
 - genotype
 - allele_counts:genotype_likelihoods:delta:MEC_variant
"""

import re
import pdb
from collections import namedtuple, OrderedDict, defaultdict

BlockInfo = namedtuple("BlockInfo", ["offset", "length", "phased",
                             "span", "MECscore", "fragments"])
Block = namedtuple("Block", ["info", "entries"])

# list of formatters - ordered dict of columns in order and functions
# that process their output
block_entry_cols = OrderedDict([("var_id", int), ("haplotype_1", int), ("haplotype_2", int),
                                ("chromosome", str), ("position", int), ("ref_allele", str), 
                                ("var_allele", str), ("genotype", lambda x: x.split(":"))])
block_entry_combined_col = OrderedDict([("allele_counts", lambda x: map(int, x.split(","))),
                                        ("genotype_likelihood", lambda x: map(float, x.split(","))),
                                          ("delta", float), ("MEC_variant", float)])
block_entry_formatters = dict(block_entry_cols.items() + block_entry_combined_col.items())

BlockEntry = namedtuple("BlockEntry", block_entry_cols.keys() + block_entry_combined_col.keys())
BLOCK_SEP = "********"
BLOCK_MATCHER = re.compile(r'BLOCK: offset: (?P<offset>\d+) len: (?P<length>\d+) phased: (?P<phased>\d+) '
                           'SPAN: (?P<span>\d+) MECscore (?P<MECscore>[\.\d]+) fragments (?P<fragments>\d+)')

last_col_matcher = ":".join(r"(?P<%s>[^:]+)" % s for s in block_entry_combined_col)
block_line_matcher = '\t'.join(r"(?P<%s>[^\t]+)" % s for s in block_entry_cols) + r"\t%s" % last_col_matcher
BLOCK_LINE_MATCHER = re.compile(block_line_matcher)

class HapCut(object):
    def __init__(self, hapcut_file):
        self.file = hapcut_file

    def blocks(self):
        """
        Iterate through phased blocks.
        """
        while True:
            line = next(self.file).strip()
            if line.startswith("BLOCK: "):
                m = BLOCK_MATCHER.match(line)
                if m is None:
                    raise ValueError("Incorrectly formatted BLOCK line.")
                block_dict = dict((k, float(v)) for k, v in m.groupdict().iteritems())
                block_info = BlockInfo(**block_dict)
                block_lines = list()
                # now, suck up other lines
                line = next(self.file).strip()
                try:
                    while not line.startswith(BLOCK_SEP):
                        m = BLOCK_LINE_MATCHER.match(line)
                        if m is None:
                            raise ValueError("Incorrectly formatted phased line.")
                        block_line_dict = dict((k, block_entry_formatters[k](v)) for k, v in m.groupdict().items())
                        this_block_line = BlockEntry(**block_line_dict)
                        block_lines.append(this_block_line)
                        line = next(self.file).strip()
                    yield Block(info=block_info, entries=block_lines)
                except StopIteration:
                    yield Block(info=block_info, entries=block_lines)
                
    def to_dict(self):
        """
        Create a dictionary of lists, each key being a reference
        sequence (chromosome or contig), and each value being a list
        of phased blocks.
        """
        self._contigs = defaultdict(list)
        for phased_block in self.blocks():
            entries = phased_block.entries
            refname = phased_block.entries[0].chromosome
            self._contigs[refname].append(phased_block)
        return self._contigs        

if __name__ == "__main__":
    a = "data/a.out"
    b = "data/hapcut-2012-08-06.out"
    h = HapCut(open(a))
    
