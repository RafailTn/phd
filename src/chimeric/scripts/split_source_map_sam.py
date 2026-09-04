#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import fire

import pysam
from seqflow import logger


@logger.catch()
def split_source_map_sam(sam, out):
    """
    Split a SAM file into 25 small SAM files.
    
    :param sam: Path to a SAM file.
    :param out: Path to an output file.
    
    Note: only SAM file > 200M will be split, but 25 small BAM files will always be created,
    even some files may not have alignments in them.
    """
    
    pysam_verbosity = pysam.set_verbosity(0)
    if os.path.getsize(sam) / 1024 / 1024 <= 200:
        pysam_verbosity = pysam.set_verbosity(0)
        with pysam.AlignmentFile(sam, 'r') as f, pysam.AlignmentFile(out, 'wb', template=f) as o:
            for read in f.fetch():
                o.write(read)
    else:
        with pysam.AlignmentFile(sam, 'r') as f:
            writers = {x + y: pysam.AlignmentFile(sam.replace('.sam', f'.{x}{y}.bam'), 'wb', template=f)
                       for x in 'ATGCN' for y in 'ATGCN'}
            for read in f.fetch():
                xy = read.query_name.split('_')[1][:2]
                writers[xy].write(read)
            _ = [v.close() for v in writers.values()]
    pysam.set_verbosity(pysam_verbosity)
    

if __name__ == '__main__':
    fire.Fire(split_source_map_sam)
