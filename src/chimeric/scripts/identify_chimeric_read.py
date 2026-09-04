#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse

import pysam
import pandas as pd
import cmder
from seqflow import logger


def _sam_to_df(sam, tag, keep_sam_string=False, bam=False, min_length=0):
    df = []
    logger.debug(f'Parsing {sam} ...')
    too_long = set()
    pysam_verbosity = pysam.set_verbosity(0)
    with pysam.AlignmentFile(sam, 'rb' if bam else 'r') as s:
        for i, read in enumerate(s.fetch(until_eof=True), start=1):
            if min_length and read.reference_length > read.query_length - min_length:
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


def _overlap(x1, x2, y1, y2):
    if x2 < y1:
        gap = y1 - x2
    elif y2 < x1:
        gap = x1 - y2
    else:
        region = list(range(x1, x2, 1))
        gap = -sum([1 if i in region else 0 for i in range(y1, y2, 1)])
    return gap


def _get_sam_header(sam):
    with open(sam) as f:
        header = ''.join([line for line in f if line.startswith('@')])
    return header


def _dedup(df, stag, tag, output=''):
    df['umi'] = df.name.str.split('_', expand=True)[1]
    logger.info(f'Before dedup: {df.shape[0]:,}')
    columns = ['umi', 'sequence', f'reference_{stag}', f'map_to_{stag}_ref_start', f'map_to_{stag}_ref_stop',
               f'reference_{tag}', f'map_to_{tag}_ref_start', f'map_to_{tag}_ref_stop']
    if output:
        dd = df[df.duplicated(subset=columns)]
        logger.info(f'Duplicated: {dd.shape[0]:,}')
        dd = dd.sort_values(by=['umi', 'sequence'])
        dd.to_csv(output, index=False, columns=columns + ['name', 'sam'])
    df = df.drop_duplicates(subset=columns)
    logger.info(f'Got {df.shape[0]:,} alignments after de-duplication.')
    df = df.drop(columns=['umi'])
    return df


@logger.catch()
def main():
    parser = argparse.ArgumentParser(description='Identify chimeric reads in BAM file.')
    parser.add_argument('bam', help='Path to a BAM file.')
    parser.add_argument('count', help='Path to count file.')
    parser.add_argument('out', help='Path to output file.')
    args = parser.parse_args()
    bam, count, out = args.bam, args.count, args.out

    _, stag, rtag, _ = bam.rsplit('.', 3)
    dd = _sam_to_df(bam, bam=True, keep_sam_string=True, tag=rtag)
    if dd.empty:
        logger.warning(f'No alignments was found in {bam}')
        return
    dd = dd.drop(columns=['read_length', 'sequence'])
    d1 = dd.groupby('name')['reference'].apply(lambda x: '|'.join(sorted(x))).reset_index()
    columns=['reference']
    d2 = dd.sort_values(by='reference').drop_duplicates(subset=['name']).drop(columns=columns)
    d3 = pd.merge(d1, d2)
    
    snorna_csv = bam.replace(f'.{rtag}', '').replace('.bam', '.target.csv').replace(stag, f'mask.map.to.{stag}')
    logger.info(f'Loading {snorna_csv} ...')
    
    df = pd.read_csv(snorna_csv)
    logger.info(f'Loaded {df.shape[0]:,} records from {snorna_csv} ...')
    
    logger.info(f'Merging {bam} [{dd.shape[0]:,}] and {snorna_csv} [{df.shape[0]:,}] ...')
    df = pd.merge(df, d3, on=['name'], how='inner', suffixes=(f'_{stag}', f'_{rtag}'))
    if df.empty:
        logger.warning(f'No candidate chimeric reads found in {bam}!')
        df['gap'] = None
    else:
        logger.info(f'Got {df.shape[0]:,} common records')
        
        logger.info(f'Calculating gaps or overlaps between partners ...')
        df['gap'] = df.apply(lambda row: _overlap(getattr(row, f'map_to_{stag}_read_start'),
                                                  getattr(row, f'map_to_{stag}_read_stop'),
                                                  getattr(row, f'map_to_{rtag}_read_start'),
                                                  getattr(row, f'map_to_{rtag}_read_stop')), axis=1)
        logger.info(f'Calculated gaps or overlaps for {df.shape[0]:,} ...')
        logger.info(f'Identifying chimeric reads in {bam} ...')
        df = df[(df['gap'].abs() <= 4) & (df[f'map_to_{rtag}_length'] >= 16)]
    
    df[f'map_to_{stag}_read_start'] = df[f'map_to_{stag}_read_start'] + 1
    df[f'map_to_{stag}_ref_start'] = df[f'map_to_{stag}_ref_start'] + 1
    df[f'map_to_{rtag}_read_start'] = df[f'map_to_{rtag}_read_start'] + 1
    df[f'map_to_{rtag}_ref_start'] = df[f'map_to_{rtag}_ref_start'] + 1
    logger.info(f'Identified {df.shape[0]:,} chimeric reads in {bam}')
    
    n1 = df['name'].unique().shape[0]
    if not df.empty:
        df = _dedup(df, stag, rtag, output=out.replace('.bam', '.duplicated.csv'))
    n2 = df['name'].unique().shape[0]
    cmder.run(f"echo {bam.replace('.bam', '')} {n1} {n2} >> {count}")
    results = {}
    start, stop = f'map_to_{stag}_read_start', f'map_to_{stag}_read_stop'
    ref_start, ref_stop = f'map_to_{stag}_ref_start', f'map_to_{stag}_ref_stop'
    target_start, target_stop = f'map_to_{rtag}_read_start', f'map_to_{rtag}_read_stop'
    target_ref_start, target_ref_stop = f'map_to_{rtag}_ref_start', f'map_to_{rtag}_ref_stop'
    df = df.rename(columns={'name': 'read_name'})
    df.to_csv(out.replace('.chimeras.bam', '.chimeras.csv'), index=False, columns=[c for c in df.columns if c != 'sam'])
    for row in df.itertuples():
        xx = {'XA': getattr(row, f'reference_{stag}'),
              'XB': f'{getattr(row, start)}_{getattr(row, stop)}_{getattr(row, ref_start)}_{getattr(row, ref_stop)}',
              'XC': f'{getattr(row, target_start)}_{getattr(row, target_stop)}_'
                    f'{getattr(row, target_ref_start)}_{getattr(row, target_ref_stop)}'}
        results[row.read_name] = xx
    pysam_verbosity = pysam.set_verbosity(0)
    with pysam.AlignmentFile(bam, 'r') as f, pysam.AlignmentFile(out, 'wb', template=f) as o:
        for read in f.fetch():
            if read.query_name in results:
                fields = results[read.query_name]
                for k, v in fields.items():
                    read.set_tag(k, v)
                o.write(read)
    pysam.set_verbosity(pysam_verbosity)
    pysam.index(out)
    

if __name__ == '__main__':
    main()
