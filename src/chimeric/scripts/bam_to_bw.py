#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse

import pysam
import cmder
from seqflow import logger


@logger.catch()
def main():
    parser = argparse.ArgumentParser(description='Make BigWig tracks using from BAM file.')
    parser.add_argument('bam', help='Path to a BAM file.')
    args = parser.parse_args()
    bam = args.bam

    reads = int(cmder.run(f'samtools view -c -F 0x4 {bam}').stdout.read())
    if reads == 0:
        logger.error(f'No mapped reads was found in {bam}, making BigWig tracks aborted.')
        return
    length = f'{bam}.to.bw.genome.length.tsv'
    cmder.run(f'samtools view -H {bam} | grep @SQ | cut -f2,3 | sed s/SN:// | sed s/LN:// > {length}')
    scale = 1000000.0 / reads
    pos_bw, neg_bw = bam.replace('.bam', '.pos.bw'), bam.replace('.bam', '.neg.bw')
    pos_bg, neg_bg = f'{pos_bw}.pos.bg', f'{neg_bw}.neg.bg'
    sort = 'LC_COLLATE=C sort -k1,1 -k2,2n'
    cmd = f'genomeCoverageBed -ibam {bam} -bg -scale {scale} -strand + -du -split | {sort} > {pos_bg}'
    cmder.run(cmd, msg=f'Calculating genome coverage for {bam} (+ strand) ...', exit_on_error=True)
    if os.path.getsize(pos_bg):
        cmd = f'bedGraphToBigWig {pos_bg} {length} {pos_bw}'
        cmder.run(cmd, msg=f'Converting {pos_bg} to {pos_bw} ...', exit_on_error=True)
    else:
        logger.debug(f'No entry was found in {pos_bg}, converting BedGraph to BigWig skipped.')
    cmder.run(f'rm {pos_bg}')
    
    cmd = f'genomeCoverageBed -ibam {bam} -bg -scale -{scale} -strand - -du -split | {sort} > {neg_bg}'
    cmder.run(cmd, msg=f'Calculating genome coverage for {bam} (- strand) ...', exit_on_error=True)
    if os.path.getsize(neg_bg):
        cmd = f'bedGraphToBigWig {neg_bg} {length} {neg_bw}'
        cmder.run(cmd, msg=f'Converting {neg_bg} to {neg_bw} ...', exit_on_error=True)
    else:
        logger.debug(f'No entry was found in {neg_bg}, converting BedGraph to BigWig skipped.')
    cmder.run(f'rm {neg_bg}')
    
    cmder.run(f'rm {length}')
    

if __name__ == '__main__':
    main()
