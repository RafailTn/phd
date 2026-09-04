#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import fire

import pysam
from seqflow import logger
from Bio import SeqIO


@logger.catch()
def compile_unmapped_putative_target(uid, tag, tags, fasta):
    """
    Compile unmapped reads into a FASTA file.
    
    :param uid: str, path to a FASTA file.
    :param tag: str, tag name for the source RNA.
    :param tags: str, comma separated tag names for target RNAs.
    :param fasta: str, path to output FASTA file.
    """
    
    mapped = set()
    for t in tags:
        sam, n = f'{uid}.{tag}.{t}.sam', 0
        logger.info(f'Loading mapped reads in {sam} ...')
        with open(sam) as f:
            for line in f:
                if not line.startswith('@') and line.strip():
                    mapped.add(line.split('\t')[0])
                    n += 1
            logger.info(f'Loaded {n:,} mapped reads from {sam}')
    
    logger.info(f'Loaded {len(mapped):,} mapped reads for sample {uid}')
    logger.info(f'Writing unmapped reads to {fasta} ...')
    reads = (read for read in SeqIO.parse(f'{uid}.mask.map.to.{tag}.true.target.fasta', 'fasta')
             if read.id not in mapped)
    n = SeqIO.write(reads, fasta, 'fasta')
    logger.info(f'Saved {n:,} reads not mapped to any RNAs.')
    

if __name__ == '__main__':
    fire.Fire(compile_unmapped_putative_target)
