#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Convert a FASTQ file to FASTA.

Replacement for the `fastq_to_fasta` console script that upstream sno-chimeras.py
imports from the VanNostrand lab's unpublished `iToolBox` package. The pipeline
calls it as two bare positional arguments::

    fastq_to_fasta INPUT.fastq[.gz] OUTPUT.fasta

which is *not* the FASTX-Toolkit interface (-i/-o), so the FASTX binary cannot be
substituted for it.

Sequences are written unwrapped on a single line, which is what every downstream
consumer in the pipeline (bowtie2 -f, Bio.SeqIO) expects.
"""

import gzip
import os
import sys

import fire
from seqflow import logger


def _opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


@logger.catch()
def fastq_to_fasta(fastq, fasta):
    """
    Convert a FASTQ file to a FASTA file.

    :param fastq: str, path to a FASTQ file, optionally gzip compressed.
    :param fasta: str, path to the output FASTA file.
    """

    if not os.path.exists(fastq):
        logger.error(f'FASTQ file {fastq} does not exist.')
        sys.exit(1)

    n = 0
    with _opener(fastq) as f, open(fasta, 'w') as o:
        for i, line in enumerate(f):
            r = i % 4
            if r == 0:
                # Drop the leading '@' and everything after the first whitespace so
                # the FASTA id matches the SAM QNAME bowtie2/STAR will emit.
                name = line[1:].split()[0] if len(line) > 1 else ''
                if not name:
                    logger.error(f'Malformed FASTQ record at line {i + 1} of {fastq}.')
                    sys.exit(1)
            elif r == 1:
                o.write(f'>{name}\n{line.strip()}\n')
                n += 1
    logger.debug(f'Converted {n:,} reads in {fastq} to FASTA file {fasta}.')


if __name__ == '__main__':
    fire.Fire(fastq_to_fasta)
