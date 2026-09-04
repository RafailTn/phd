#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import fire

import pysam
from seqflow import logger


@logger.catch()
def sam_to_bam(sam, bam, name_reset):
    """
    Convert a SAM file to BAM and reset read name if needed.
    
    :param sam: str, path to a SAM file.
    :param bam: str, path to the output BAM file.
    :param name_reset: bool, True indicates the read name needs to be reset,
        False will ignore read name reset
    """

    tmp_bam = sam.replace('.sam', '.tmp.bam')
    pysam_verbosity = pysam.set_verbosity(0)
    with pysam.AlignmentFile(sam) as s, pysam.AlignmentFile(tmp_bam, "wb", template=s) as o:
        for read in s.fetch():
            if name_reset:
                name, offset = read.query_name.rsplit('_', maxsplit=1)
                read.set_tag('YN', int(offset))
                read.query_name = name
            o.write(read)
    pysam.set_verbosity(pysam_verbosity)
    pysam.sort(tmp_bam, "-o", bam)
    pysam.index(bam)
    os.unlink(tmp_bam)
    

if __name__ == '__main__':
    fire.Fire(sam_to_bam)
