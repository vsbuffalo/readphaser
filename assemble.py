import sys
import re
import pdb
from time import sleep
from multiprocessing import Queue, Pool, Process
from Queue import Empty
import argparse
import fermi as fm
import gzip
from collections import namedtuple
from readfq import readfq

PhasedBlock = namedtuple("PhasedBlock", ["contig", "block", "phase"])
FASTQEntry = namedtuple("FASTQEntry", ["name", "seq", "qual"])

def contig_writer(queue_out, file_handle):
    while True:
        try:
            tigs = queue_out.get(block=False)
            if tigs is None:
                break
            for tig in tigs:
                file_handle.write("@%s\n%s\n+\n%s\n" % (tig.header, tig.seq, tig.qual))
        except Exception, e:
            if isinstance(e, Empty):
                continue
            sys.stderr.write("[error] exception in queue_writer: %s" % e)
            break


def assembly(block, readlist):
    # fermi = fm.Fermi()
    # for read in readlist:
    #     fermi.addseq(read.seq, read.qual)
    # fermi.correct()
    # tigs = fermi.assemble(unitig_k=-1, do_clean=True)
    # if tigs is None:
    #     #sys.stdout.write("%s\tNA")
    #     return None
    # root_name = "%s %s" % (block.contig[2:], " ".join((block.block, block.phase)))
    # tigs = fermi.fastq_to_list(tigs, root_name)
    # return tigs
    sleep(3)
    return [FASTQEntry("test", "ATGC", "BBBB")]

def assembly_runner(queue_in, queue_out):
    while True:
        try:
            readlist = queue_in.get(block=False)
            if readlist is None:
                queue_out.put(None)
                break
            block = readlist.pop(0)
            sys.stderr.write("[assembly_runner] processing %s\n" % " ".join((block.contig, block.block, block.phase)))
            tigs = assembly(block, readlist)
            if tigs is not None:
                queue_out.put(tigs)
        except Exception, e:
            if isinstance(e, Empty):
                continue
            pdb.set_trace()
            sys.stderr.write("[error] exception in assembly_runner: %s" % e)
            break

def parse_header(header):
    readname, contig, block, phase = header.strip().split()
    return readname, PhasedBlock(contig, block, phase)
    
def FASTQ_to_readlists(file_handle):
    past_blocks = set()
    last_contig, last_block, last_phase = None, None, None
    # we must assert that everything is in order TODO
    readlist = list()
    for header, seq, qual in readfq(file_handle):
        name, block = parse_header(header)
        block_key = "-".join((block.contig, block.block, block.phase))
        same = (last_contig == block.contig, last_block == block.block, last_phase == block.phase)
        is_first = None in (last_contig, last_block, last_phase)
        if not is_first and not all(same):
            assert(block_key not in past_blocks)
            if len(readlist) > 1: # only block
                yield readlist
                            
            #sys.stdout.write("\t") # TODO stats

            past_blocks.add(block_key)
            readlist = list()
            # first item of every readlist is block info
            readlist.append(block)
        if is_first:
            readlist.append(block)
        readlist.append(FASTQEntry(header, seq, qual))
        last_contig, last_block, last_phase = block.contig, block.block, block.phase
    yield readlist

def assembly_queuer(queue, file_handle):
    for readlist in FASTQ_to_readlists(file_handle):
        sys.stderr.write("[assembly_queuer] adding %s to queue\n" % ' '.join(readlist[0]))
        queue.put(readlist)
    queue.put(None)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.stderr.write("usage: assemble.py phased-reads.fastq num_processes output_filename.fastq")
        sys.exit(1)

    filename = sys.argv[1]
    openfun = gzip.GzipFile if filename.endswith(".gz") else open
    num_processes = int(sys.argv[2])
    input_filehandle = openfun(filename, 'r')
    output_filehandle = open(sys.argv[3], 'w')

    assembly_queue = Queue()
    writer_queue = Queue()

    if num_processes == 1:
        assembly_queuer(assembly_queue, input_filehandle)
        assembly_runner(assembly_queue, writer_queue)
        contig_writer(writer_queue, output_filehandle)
    else:    
        assembly_queuer_processes = Process(target=assembly_queuer, args=(assembly_queue, input_filehandle))
        assembly_processes = [Process(target=assembly_runner, args=(assembly_queue, writer_queue)) for _ in range(num_processes)]
        writer_processes = Process(target=contig_writer, args=(writer_queue, output_filehandle))

        sys.stderr.write("[main] starting writer_process\n")
        writer_processes.start()

        sys.stderr.write("[main] starting assembly_queuer\n")
        assembly_queuer_processes.start()
        sys.stderr.write("[main] joining assembly_queuer\n")
        assembly_queuer_processes.join()

        for i, p in enumerate(assembly_processes):
            sys.stderr.write("[main] starting assembly_process %d\n" % i)
            p.start()
        for i, p in enumerate(assembly_processes):
            sys.stderr.write("[main] joining assembly_process %d\n" % i)
            p.join()

        sys.stderr.write("[main] joining writer_process\n")
        writer_processes.join()
