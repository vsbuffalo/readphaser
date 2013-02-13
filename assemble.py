import itertools
import sys
from datetime import datetime
import re
import pdb
from time import sleep
from multiprocessing import Queue, Pool, Process, Manager
from Queue import Empty
import argparse
import fermi as fm
import gzip
from collections import namedtuple
from readfq import readfq

PhasedBlock = namedtuple("PhasedBlock", ["contig", "block", "phase"])
FASTQEntry = namedtuple("FASTQEntry", ["name", "seq", "qual"])

# def contig_writer(queue_out, file_handle):
#     while True:
#         try:
#             tigs = queue_out.get(block=False)
#             if tigs is None:
#                 break
#             for tig in tigs:
#                 file_handle.write("@%s\n%s\n+\n%s\n" % (tig.header, tig.seq, tig.qual))
#         except Exception, e:
#             if isinstance(e, Empty):
#                 continue
#             sys.stderr.write("[error] exception in queue_writer: %s" % e)
#             break


def assembly(readlist):
    start = datetime.now().isoformat()
    block = readlist[0]
    fermi = fm.Fermi()
    for read in readlist[1:]:
        name, seq, qual = read
        fermi.addseq(seq, qual)
    
#    pdb.set_trace()
    correct_out = fermi.correct()
    tigs = fermi.assemble(unitig_k=int(-1), do_clean=True)
    if tigs is None:
        #sys.stdout.write("%s\tNA")
        return None
    root_name = "-".join((block.contig[3:], block.block, block.phase))
    tigs = fermi.fastq_to_list(tigs, root_name)
    end = datetime.now().isoformat()
    sys.stdout.write("[assembly] completed for %s (s: %s e: %s)\n" % (" ".join(block), start, end))
    return tigs

# def assembly_runner(queue_in, queue_out):
#     while True:
#         try:
#             readlist = queue_in.get(block=False)
#             if readlist is None:
#                 queue_out.put(None)
#                 break
#             block = readlist.pop(0)
#             sys.stderr.write("[assembly_runner] processing %s\n" % " ".join((block.contig, block.block, block.phase)))
#             tigs = assembly(block, readlist)
#             if tigs is not None:
#                 queue_out.put(tigs)
#         except Exception, e:
#             if isinstance(e, Empty):
#                 continue
#             pdb.set_trace()
#             sys.stderr.write("[error] exception in assembly_runner: %s" % e)
#             break

def parse_header(header):
    readname, block, phase, contig = header.strip().split()
    return readname, PhasedBlock(contig, block, phase)
    
def FASTQ_to_readlists(file_handle):
    past_blocks = set()
    last_contig, last_block, last_phase = None, None, None
    # we must assert that everything is in order TODO
    readlist = list()
    for header, seq, qual in readfq(file_handle):
        name, block = parse_header(header)
        is_first = None in (last_contig, last_block, last_phase)
        if is_first:
            readlist.append(block)

        block_key = "-".join((block.contig, block.block, block.phase))
        same = (last_contig == block.contig, last_block == block.block, last_phase == block.phase)
        if not is_first and not all(same):
            assert(block_key not in past_blocks)
            if len(readlist) > 1: # not only block
                yield readlist
                            
            #sys.stdout.write("\t") # TODO stats

            past_blocks.add(block_key)
            readlist = list()
            # first item of every readlist is block info
            readlist.append(block)
        readlist.append(FASTQEntry(header, seq, qual))
        last_contig, last_block, last_phase = block.contig, block.block, block.phase
    yield readlist

# def assembly_queuer(queue, file_handle):
#     for readlist in FASTQ_to_readlists(file_handle):
#         sys.stderr.write("[assembly_queuer] adding %s to queue\n" % ' '.join(readlist[0]))
#         queue.put(readlist)
#     queue.put(None)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.stderr.write("usage: assemble.py phased-reads.fastq num_processes output_filename.fastq")
        sys.exit(1)

    filename = sys.argv[1]
    openfun = gzip.GzipFile if filename.endswith(".gz") else open
    num_processes = int(sys.argv[2])
    input_filehandle = openfun(filename, 'r')
    output_filehandle = open(sys.argv[3], 'w')

    if num_processes == 1:
        contigs_out = list()
    items = list(FASTQ_to_readlists(input_filehandle))
    for readlist in items:
        contigs_out.append(assembly(readlist))
    # else:    
    #     pool = Pool(num_processes)
    #     #contigs_out = pool.map(assembly, FASTQ_to_readlists(input_filehandle))
    #     readlists_gen = FASTQ_to_readlists(input_filehandle)
    #     contigs_out = list()
    #     contigs_out.extend(pool.imap(assembly, itertools.islice(readlists_gen, None)))
                

    for contigs in contigs_out:
        if contigs is None:
            continue
        for tig in contigs:
            output_filehandle.write("@%s\n%s\n+\n%s\n" % (tig.header, tig.seq, tig.qual))
        
        
    # readlists_gen = FASTQ_to_readlists(input_filehandle)
    # for readslist in readlists_gen:
    #     for block in readslist[1:]:
    #         print "@%s\n%s\n+\n%s" % (block.name, block.seq, block.qual)
