#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import fire

import pysam
import pandas as pd
from seqflow import logger


def sam_to_df(sam, tag, keep_sam_string=False, bam=False, min_length=0):
    df = []
    logger.debug(f'Parsing {sam} ...')
    too_long = set()
    pysam_verbosity = pysam.set_verbosity(0)
    with pysam.AlignmentFile(sam, 'rb' if bam else 'r') as s:
        for i, read in enumerate(s.fetch(until_eof=True), start=1):
            if min_length and min_length > read.query_length - read.reference_length:
                too_long.add(read.query_name)
                continue
            offset = read.get_tag('YN') if read.has_tag('YN') else 0
            # The tag appended by STAR need to be removed
            alignment = {
                'name': read.query_name.split('/')[0],
                'sequence': read.query_sequence,
                'reference': read.reference_name,
                'read_length': read.query_length,
                f'map_to_{tag}_length': read.reference_length,
                f'map_to_{tag}_read_start': read.query_alignment_start + offset,
                f'map_to_{tag}_read_stop': read.query_alignment_end + offset,
                f'map_to_{tag}_ref_start': read.reference_start,
                f'map_to_{tag}_ref_stop': read.reference_end,
                f'map_to_{tag}_score': read.get_tag('AS') if read.has_tag('AS') else 0,
                f'map_to_{tag}_strand': '-' if read.is_reverse else '+'
            }
            if keep_sam_string:
                alignment['sam'] = read.to_string()
            df.append(alignment)
    pysam.set_verbosity(pysam_verbosity)
    if too_long:
        logger.debug(f'Found {len(too_long):,} reads have too many bases mapped!')
        df = [row for row in df if row['name'] not in too_long]
    df = pd.DataFrame(df)
    logger.debug(f'Parsing {sam} complete: {df.shape[0]:,} alignments were loaded into DataFrame.')
    return df


def _populate_target(row, tag=''):
    left, right = row['sequence'][:row[f'map_to_{tag}_read_start']], row['sequence'][row[f'map_to_{tag}_read_stop']:]
    target_seq, offset = (left, 0) if len(left) >= len(right) else (right, row[f'map_to_{tag}_read_stop'])
    return (f'{row["name"]}_{offset}', target_seq) if len(target_seq) >= 16 else ('', '')


@logger.catch()
def find_putative_target(bam, fasta, tag, debug):
    """
    Find putative target using a BAM file and a FASTA file.
    
    :param bam: str, path to a BAM file.
    :param fasta: str, path to a FASTA file.
    :param tag: str, tag name for the source RNA.
    :param debug: bool, True or False for whether to invoke debug mode to save intermedia files.
    """

    if os.path.exists(bam):
        logger.debug(f'Loading alignments from {bam} ...')
        # Both aligned and unaligned length were set to 18 and target min length was set to 16.
        # This will allow source and target have at least 2-base gap
        df = sam_to_df(bam, tag, bam=True, min_length=18)
        if df.empty:
            with open(fasta, 'w') as o:
                o.write('')
            logger.debug('No alignments were found, empty target fasta file was created.')
            return
        
        logger.debug(f'Filtering alignments by read name and mapping score ...')
        df = df[df.groupby(by='name')[f'map_to_{tag}_score'].transform('max') == df[f'map_to_{tag}_score']]
        logger.debug(f'Got {df.shape[0]:,} alignments after filtered by name and score.')
        if df.empty:
            return
        
        logger.debug(f'Filtering alignments by read name, aligned reference, and aligned length ...')
        length = f'map_to_{tag}_length'
        df = df[df.groupby(by=['name', 'reference'])[length].transform('max') == df[length]]
        logger.debug(f'Got {df.shape[0]:,} alignments after filtered by name, reference, and aligned length.')
        if df.empty:
            return
        
        logger.info(f'Filtering alignments by aligned length ...')
        df = df[df[length] >= 18]
        logger.debug(f'Got {df.shape[0]:,} alignments with valid length for retrieving putative targets.')
        if df.empty:
            return
        
        logger.debug('Determining putative targets ...')
        df[['target_name', 'target_sequence']] = df.apply(_populate_target, axis=1, tag=tag, result_type='expand')
        df = df[df.target_name != '']
        logger.debug(f'Got {df.shape[0]:,} alignments with putative targets.')
        if df.empty:
            return
        
        d1 = df.groupby('name')['reference'].apply(lambda x: '|'.join(sorted(x))).reset_index()
        columns = ['reference', 'read_length', f'map_to_{tag}_length', f'map_to_{tag}_score', f'map_to_{tag}_strand']
        d2 = df.sort_values(by='reference').drop_duplicates(subset=['name']).drop(columns=columns)
        df = pd.merge(d1, d2)
        
        logger.debug('Saving putative targets to file ...')
        with open(fasta, 'w') as o:
            o.writelines(f'>{row.target_name}\n{row.target_sequence}\n' for row in df.itertuples())
        logger.debug(f'Saved {df.shape[0]:,} distinct putative target sequences to {fasta}.')
        
        df = df.drop(columns=['target_name', 'target_sequence'])
        df.to_csv(fasta.replace('.fasta', '.csv'), index=False)
    else:
        with open(fasta, 'w') as o:
            o.write('')
        logger.debug(f'BAM file {bam} was not found, empty target fasta file was created.')
    

if __name__ == '__main__':
    fire.Fire(find_putative_target)
