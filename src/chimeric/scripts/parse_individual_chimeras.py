#!/usr/bin/env python

"""
Separate chimeras for RNAs meet the minimum number of chimeras into individual BAM files.
"""
import argparse
import glob
import sys
from collections import defaultdict

import pysam
import cmder
from seqflow import logger


@logger.catch()
def main():
    parser = argparse.ArgumentParser(prog='split_individual_chimeras', description=__doc__)
    parser.add_argument('bam', help='Path to a chimeric BAM file.', type=cmder.filename)
    parser.add_argument('--uid', help='Unique identifier the chimeric BAM file associates with.')
    parser.add_argument('--tag', help='Target tag the chimeric BAM file associates with.')
    parser.add_argument('--species', help='Species passes to clipper as species for peak calling, '
                                          'without specifying, no peaks will be called.')
    parser.add_argument('--minimum_chimeras',
                        help='Individual snoRNA with at least this minimum number of chimeric reads '
                             'will be separated into individual files, default: %(default)s.', default=1000, type=int)
    parser.add_argument('--no_bw', action='store_true', help='Flag that disables creation of BigWig files.')
    parser.add_argument('--cpus', type=int, default=8,
                        help='Number of CPU cores can be used for paralleling your job, default: %(default)s.')
    parser.add_argument('--count_output', help='Output file for saving counts of individual chimeras.')
    
    args = parser.parse_args()
    if not args.count_output:
        if bam.endswith('.chimeras.bam'):
            setattr(args, 'count_output', arg.bam.replace('.chimeras.bam', '.individual.chimeras.count.csv'))
        else:
            setattr(args, 'count_output', f'{bam}.individual.chimeras.count.csv')
    setattr(args, 'uid', args.uid or args.bam)
    setattr(args, 'tag', args.tag or args.bam)

    chimeras, counts = defaultdict(list), []
    with pysam.AlignmentFile(args.bam, 'rb') as f:
        for read in f.fetch():
            source = read.get_tag('XA').split('|')
            for s in source:
                chimeras[s].append(read)
    
        if chimeras:
            for k, v in chimeras.items():
                if len(v) >= args.minimum_chimeras:
                    counts.append([args.uid, k, args.tag, str(len(v))])
                    if args.bam.endswith('.chimeras.bam'):
                        out = args.bam.replace('.chimeras.bam', f'.{k}.chimeras.bam')
                        bed = args.bam.replace('.chimeras.bam', f'.{k}.chimeras.peak.bed')
                    else:
                        out = f'{args.bam}.{k}.chimeras.bam'
                        bed = f'{args.bam}.{k}.chimeras.peak.bed'
                    logger.info(f'Saving {len(v)} chimeras to {out}.')
                    with pysam.AlignmentFile(out, 'wb', template=f) as o:
                        for x in v:
                            o.write(x)
                    pysam.index(out)
                    if args.species:
                        cmd = f'clipper --species {args.species} --processors {args.cpus} --bam {out} --outfile {bed}'
                        cmder.run(cmd, msg=f'Calling peaks in {out} ...', exit_on_error=True)
                        
                    if args.no_bw:
                        pass
                    else:
                        cmder.run(f'bam_to_bw {out}', msg=f'Creating BigWig tracks using {out} ...', exit_on_error=True)
            else:
                logger.warning(f'No individual chimeras in {args.bam} meet the minimum chimeras count.')

    if counts:
        with open(args.count_output, 'a') as o:
            o.write('uID,snoRNA,target,chimeras_count\n')
            o.writelines(f'{",".join(count)}\n' for count in counts)
        
        
if __name__ == '__main__':
    main()
