#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import shutil
import subprocess
import sys
import gzip
from collections import defaultdict, Counter

import cmder
from seqflow import task, logger, Flow

# Reference directory holding the adapter FASTAs, overridable so the same
# checkout runs elsewhere without editing this file.
REF_DIR = os.environ.get('CHIMERIC_REF_DIR',
                         os.path.join(os.path.dirname(os.path.dirname(
                             os.path.dirname(os.path.abspath(__file__)))), 'ref', 'chimeric'))


parser = argparse.ArgumentParser(description='Identify snoRNA chimeric reads using bowtie and STAR.')
parser.add_argument('fastq', metavar='FASTQ', type=cmder.filename, nargs='+',
                    help='Path to one or multiple FASTQ files (multiple files separated by space)')
parser.add_argument('--uid', nargs='*',
                    help='Unique identifier for each FASTQ file(s), without specifying, basename of each FASTQ '
                         'file without the .fastq.gz extension will be used.')
parser.add_argument('--strategy', default='chimeric', choices=('eCLIP', 'modified_eCLIP', 'chimeric'),
                    help="Adapter trimming strategy, optional."
                         " eCLIP   - Two rounds trimming with overlap set to 1 and 5 for fist and second round."
                         " modified_eCLIP - Three rounds trimming with overlap set to 1 and 5 for fist and "
                         "second round, and 10 bases trimmed from 3' end for third round."
                         " chimeric - Two round trimming with overlap set to 1 for 3 times first round trimming, "
                         "and 10 bases trimmed from 3' end for second round.")
parser.add_argument('--adapters_fasta', type=cmder.filename,
                    help="Path to the FASTA file contains adapters and their sequences, "
                         "default: auto populate based on trimming strategy.")
parser.add_argument('--umi_pattern', help="Pattern of UMI for umi-tools using it to extract UMIs, "
                                          "%(default)s.", default='NNNNNNNNNN')
parser.add_argument('--min_length', help="Minimum length of reads after adapter removal, default: %(default)s.",
                    type=int, default=24)
parser.add_argument('--randomer_length', help="Length of 3 premier randomer, default: %(default)s.",
                    type=int, default=10)
parser.add_argument('--species', help='Short name for species, e.g., hg19, mm10, ..., default: %(default)s.',
                    choices=('hg19', 'mm10', 'hg19chr19', 'hg38'), default='hg19')
parser.add_argument('--source_rna_fasta', required=True, type=cmder.filename,
                    help='Path to a FASTA file contains sequences for source RNA, required.')
parser.add_argument('--source_rna_tag', help='Name or tag for source RNA, required.', required=True)
parser.add_argument('--target_rna_fasta', nargs='+', type=cmder.filename, required=True,
                    help='Path(s) to FASTA file(s) that contains sequences for target RNA(s), required.')
parser.add_argument('--target_rna_tag', nargs='+', help='Name(s) or tag(s) for target RNA(s), required.', required=True)
parser.add_argument('--source_rna_bed', help='Bed file of source RNA for filtering out false genomic chimeras.',
                    type=cmder.filename)
parser.add_argument('--target_rRNA_tag', help='Name or tag for target rRNA,default: %(default)s.', default='rRNA')
parser.add_argument('--clipper_rna_species', help='Species passes to clipper as species for rRNA chimeras '
                                                  'peak calling, otherwise, no rRNA chimeras peaks will be called.')
parser.add_argument('--minimum_chimeras', help='Individual snoRNA with at least this minimum number of chimeric '
                                               'reads will be separated into individual files, default: %(default)s.',
                    default=1000, type=int)
parser.add_argument('--genome_star_index', type=cmder.dirname,
                    help='Path to reference genome STAR index, default: auto populate based on species.')
parser.add_argument('--repeat_star_index', type=cmder.dirname,
                    help='Path to repeat element STAT index, default: auto populate based on species.')
parser.add_argument('--outdir',
                    help="Path to the output directory. Default to the current work "
                         "directory and if the specified path does not exist, it will "
                         "try to create it first.")
parser.add_argument('--cpus', type=int, default=8, 
                    help='Number of CPU cores can be used for paralleling your job, default: %(default)s.')
parser.add_argument('--outFilterMismatchNmax', type=int, default=10,
                    help='Reads only be used if it has no more mismatches than this value, default: %(default)s.')
parser.add_argument('--debug', action='store_true',
                    help='Invoke debug mode and keep all intermedia files along with the result files.')
parser.add_argument('--mask_repeat_genome_read', action='store_true',
                    help='Invoke the task that mask reads that mapped to repeat elements and genomic regions '
                         'before mapping reads to source snoNRA.')
parser.add_argument('--keep', action='store_true', help='Flag that will keep all intermedia files undeleted.')
parser.add_argument('--dryrun', action='store_true',
                    help='Print out steps and files involved in each step without actually '
                         'running the pipeline.')

args = parser.parse_args()
setattr(args, 'outdir', args.outdir or os.getcwd())
try:
    os.makedirs(args.outdir, exist_ok=True)
except OSError as e:
    logger.error(f'Failed to create outdir, {e}.')
    sys.exit(1)
os.chdir(args.outdir)

if args.uid:
    if len(args.uid) == len(args.fastq):
        FASTQ2UID = {fastq: name for fastq, name in zip(args.fastq, args.uid)}
    else:
        logger.error('The number of FASTQ files does not math the number of names.')
        sys.exit(1)
else:
    FASTQ2UID = {fastq: os.path.basename(fastq).replace('.fastq.gz', '') for fastq in args.fastq}
IDS = list(FASTQ2UID.values())

if len(args.target_rna_fasta) == len(args.target_rna_tag):
    TARGETS = {tag: fasta for tag, fasta in zip(args.target_rna_tag, args.target_rna_fasta)}
else:
    logger.error('The number of target RNA FASTA files does not equal to the number of target RNA tags.')
    sys.exit(1)

# Upstream defaulted these to indices on the van Nostrand lab cluster. Those
# paths exist on no machine here, so a missing --genome_star_index used to fail
# with someone else's absolute path; require the flags instead.
GENOMES, REPEATS = {}, {}
if not args.genome_star_index:
    logger.error('No STAR genome index: pass --genome_star_index '
                 '(run_chimeras.sh does, from --genome-index).')
    sys.exit(1)
if not args.repeat_star_index:
    logger.error('No STAR repeat index: pass --repeat_star_index '
                 '(run_chimeras.sh does, from --repeat-index).')
    sys.exit(1)
GENOMES[args.species] = args.genome_star_index
# CLIPPER_RNA_SPECIES = {
#     'hg19': 'hg19.rRNA.45S',
#     'hg19chr19': 'hg19.rRNA.45S',
#     'hg38': 'hg19.rRNA.45S',
#     'mm10': 'mm10.rRNA.5S.45S'
# }
# setattr(args, 'clipper_rna_species', args.clipper_rna_species or CLIPPER_RNA_SPECIES[args.species])
REPEATS[args.species] = args.repeat_star_index
if args.target_rRNA_tag not in args.target_rna_tag:
    logger.error(f'Target rRNA tag "{args.target_rRNA_tag}" not found in target rna tag {args.target_rna_tag}.')
    sys.exit(1)

STAG, TAG, GTAG = args.source_rna_tag, args.target_rna_tag, args.species
TAGS = TAG + [args.species]
COUNT = f'{STAG}.{".".join(TAG)}.chimeras.count.tsv'

if not args.dryrun:
    script = f'{".".join(IDS)}.snorna.chimeras.pipeline.sh'
    if os.path.isfile(script):
        os.unlink(script)
    cmder.logger.add(script, format="{message}",
                     filter=lambda record: record["level"].name == "DEBUG" and not record['message'].startswith('Task'))
    cmder.logger.debug(f'#!/usr/bin/env bash')
    cmder.logger.debug(f'set -e')
    cmder.logger.debug(f'[ -d {args.outdir} ] || mkdir {args.outdir}')
    cmder.logger.debug(f'cd {args.outdir}')


@task(inputs=args.fastq, outputs=lambda i: f'{FASTQ2UID[i]}.trim.fastq.gz')
def extract_umi_and_cut_adapter(fastq, output):
    """Extract UMIs and trim off adapters."""

    if output.endswith('.trim.fastq.gz'):
        log = output.replace('.trim.fastq.gz', '.cut.adapt.log')
    else:
        log = f'{output}.cut.adapt.log'
    if args.strategy == 'eCLIP':
        adapters_fasta = os.path.join(REF_DIR, 'a_adapters.fasta')
        cmd = ['umi_tools', 'extract', '--random-seed', 1, '--bc-pattern', args.umi_pattern,
               '--stdin', fastq, '--log', '/dev/null', '|',

               'cutadapt', '-j', args.cpus, '--match-read-wildcards',
               '--times', 1, '-e', 0.1, '--quality-cutoff', 6, '-m', args.min_length,
               '-a', f'file:{args.adapters_fasta or adapters_fasta}', '-O', 1, '-',
               '--report', 'minimal', '2>', log, '|',

               'cutadapt', '-j', args.cpus, '--match-read-wildcards',
               '--times', 1, '-e', 0.1, '--quality-cutoff', 6, '-m', args.min_length,
               '-a', f'file:{args.adapters_fasta or adapters_fasta}', '-O', 5, '-', '-o', output,
               '--report', 'minimal', '>>', log, '2>', '/dev/null']
    elif args.strategy == 'modified_eCLIP':
        adapters_fasta = os.path.join(REF_DIR, 'se.adapters.fasta')
        cmd = ['umi_tools', 'extract', '--random-seed', 1, '--bc-pattern', args.umi_pattern,
               '--stdin', fastq, '--log', '/dev/null', '|',

               'cutadapt', '-j', args.cpus, '--match-read-wildcards',
               '--times', 1, '-e', 0.1, '--quality-cutoff', 6, '-m', args.min_length + args.randomer_length,
               '-a', f'file:{args.adapters_fasta or adapters_fasta}', '-O', 1, '-',
               '--report', 'minimal', '2>', log, '|',

               'cutadapt', '-j', args.cpus, '--match-read-wildcards',
               '--times', 1, '-e', 0.1, '--quality-cutoff', 6, '-m', args.min_length + args.randomer_length,
               '-a', f'file:{args.adapters_fasta or adapters_fasta}', '-O', 5, '-',
               '--report', 'minimal', '2>>', log, '|',

               'cutadapt', '-j', args.cpus, '-u', -1 * args.randomer_length, '-', '-o', output,
               '--report', 'minimal', '>>', log, '2>', '/dev/null']
    elif args.strategy == 'chimeric':
        adapters_fasta = os.path.join(REF_DIR, 'se.2.round.adapters.fasta')
        cmd = ['umi_tools', 'extract', '--random-seed', 1, '--bc-pattern', args.umi_pattern,
               '--stdin', fastq, '--log', '/dev/null', '|',

               'cutadapt', '-j', args.cpus, '--match-read-wildcards',
               '--times', 3, '-e', 0.1, '--quality-cutoff', 6, '-m', args.min_length + args.randomer_length,
               '-a', f'file:{args.adapters_fasta or adapters_fasta}', '-O', 1, '-',
               '--report', 'minimal', '2>', log, '|',

               'cutadapt', '-j', args.cpus, '-u', -1 * args.randomer_length, '-o', output, '-',
               '--report', 'minimal', '>>', log, '2>', '/dev/null']
    else:
        logger.error(f'Invalid cut adapter strategy: {args.strategy}.')
        sys.exit(1)

    p = cmder.run(cmd, msg=f'Extracting UMIs and trimming adapters from {fastq} ...',
                  pmt=args.debug, debug=args.debug, exit_on_error=True)
    if p.returncode:
        sys.exit(p.returncode)


@task(inputs=extract_umi_and_cut_adapter, outputs=lambda i: i.replace('.trim.fastq.gz', '.mask.fasta'))
def mask_repeat_genome_read(fastq, out):
    if args.mask_repeat_genome_read:
        prefix = fastq.replace('.trim.fastq.gz', '.repeat.map')
        mate = fastq.replace('.trim.fastq.gz', '.repeat.unmap.fastq')
        n = sum([int(os.path.exists(f'{name}.repeat.unmap.fastq')) for name in IDS])
        genome_load = 'LoadAndRemove' if n == len(IDS) - 1 else 'LoadAndKeep'
        try:
            os.makedirs(prefix, exist_ok=True)
            cmder.logger.debug(f'[ -d {prefix} ] || mkdir {prefix}')
            cmd = ['STAR',
                   '--runMode', 'alignReads',
                   '--runThreadN', args.cpus,
                   '--alignEndsType', 'EndToEnd',
                   '--genomeDir', REPEATS[args.species],
                   '--genomeLoad', genome_load,
                   '--outBAMcompression', 10,
                   '--outFileNamePrefix', f"{prefix}/",
                   '--outFilterMultimapNmax', 100,
                   '--outFilterMultimapScoreRange', 1,
                   '--outFilterScoreMin', 10,
                   '--outFilterMismatchNmax', args.outFilterMismatchNmax,
                   '--outFilterType', 'BySJout',
                   '--outReadsUnmapped', 'Fastx',
                   '--outSAMattrRGline', 'ID:foo',
                   '--outSAMattributes', 'All',
                   '--outSAMmode', 'Full',
                   '--outSAMtype', 'BAM', 'Unsorted',
                   '--outSAMunmapped', 'None',
                   '--outStd', 'Log',
                   '--readFilesCommand', 'zcat',
                   '--readFilesIn', fastq]

            cmder.run(cmd, msg=f'Map reads in {fastq} to repeat elements ...',
                      pmt=args.debug, debug=args.debug, exit_on_error=True)
            if args.debug:
                cmder.run(f'cp {prefix}/Log.final.out {mate.replace(".repeat.unmap.fastq", ".mask.repeat.map.log")}')
                cmder.run(f'cp {prefix}/Unmapped.out.mate1 {mate}')
            else:
                cmder.run(f'mv {prefix}/Log.final.out {mate.replace(".repeat.unmap.fastq", ".mask.repeat.map.log")}')
                cmder.run(f'mv {prefix}/Unmapped.out.mate1 {mate}')
        finally:
            if args.keep or args.debug:
                pass
            else:
                cmder.run(f'rm -r {prefix}')

        prefix = fastq.replace('.trim.fastq.gz', '.genome.map')
        try:
            os.makedirs(prefix, exist_ok=True)
            cmder.logger.debug(f'[ -d {prefix} ] || mkdir {prefix}')
            cmd = ['STAR',
                   '--runMode', 'alignReads',
                   '--runThreadN', args.cpus,
                   '--alignEndsType', 'EndToEnd',
                   '--genomeDir', GENOMES[args.species],
                   '--genomeLoad', genome_load,
                   '--outBAMcompression', 10,
                   '--outFileNamePrefix', f"{prefix}/",
                   '--outFilterMultimapNmax', 100,
                   '--outFilterMultimapScoreRange', 1,
                   '--outFilterScoreMin', 10,
                   '--outFilterMismatchNmax', args.outFilterMismatchNmax,
                   '--outFilterType', 'BySJout',
                   '--outReadsUnmapped', 'Fastx',
                   '--outSAMattrRGline', 'ID:foo',
                   '--outSAMattributes', 'All',
                   '--outSAMmode', 'Full',
                   '--outSAMtype', 'BAM', 'Unsorted',
                   '--outSAMunmapped', 'None',
                   '--outStd', 'Log',
                   '--readFilesIn', mate]

            cmder.run(cmd, msg=f'Map repeat elements unmapped reads in {mate} to reference genome ...',
                      pmt=args.debug, debug=args.debug, exit_on_error=True)
            exe = 'cp' if args.debug else 'mv'
            cmder.run(f'{exe} {prefix}/Log.final.out {mate.replace(".repeat.unmap.fastq", ".mask.genome.map.log")}')
            cmder.run(f'fastq_to_fasta {prefix}/Unmapped.out.mate1 {out}')
        finally:
            if args.keep or args.debug:
                pass
            else:
                cmder.run(f'rm -r {prefix}')
    else:
        cmder.run(f'fastq_to_fasta {fastq} {out}', msg=f'Convert {fastq} to {fastq} '
                                                       f'(since task mask_repeat_genome_read was not enabled.')
    return out


def bowtie_index(basename, fasta):
    if not os.path.exists(f'{basename}.rev.2.ebwt'):
        cmder.run(f'bowtie-build --quiet --offrate 2 {fasta} {basename} 2> /dev/null',
                  msg=f'Building bowtie index using {fasta} ...', fmt_cmd=False,
                  pmt=args.debug, debug=args.debug, exit_on_error=True)
    return basename


def bowtie2_index(basename, fasta):
    if not os.path.exists(f'{basename}.rev.2.bt2'):
        cmder.run(f'bowtie2-build --quiet --offrate 2 {fasta} {basename} 2> /dev/null',
                  msg=f'Building bowtie2 index using {fasta} ...', fmt_cmd=False,
                  pmt=args.debug, debug=args.debug, exit_on_error=True)
    return basename
    

@task(inputs=mask_repeat_genome_read,
      outputs=lambda i: i.replace(".mask.fasta", f".mask.map.to.{STAG}.sam"))
def map_mask_read_to_source_rna(fasta, sam):
    db = bowtie2_index(f'{STAG}.bowtie2.index', args.source_rna_fasta)
    cmd = ['bowtie2', '-D', '20', '-R', '3', '-N', '0', '-L', '16', '--local', '--norc', '-i', 'S,1,0.50',
           '--score-min', 'L,16,0', '--ma', '1', '--np', '0', '--mp', '2,2', '--rdg', '5,1', '--rfg', '5,1',
           '-p', args.cpus, '--no-unal', '-x', db, '-f', '-U', fasta, '-a', '-S', sam,
           '2>', sam.replace('.sam', '.log')]
    cmder.run(cmd, msg=f'Mapping {os.path.basename(fasta)} to {STAG} ...', 
              pmt=args.debug, debug=args.debug, exit_on_error=True)


@task(inputs=map_mask_read_to_source_rna, outputs=lambda i: i.replace('.sam', '.AA.bam'), cpus=len(FASTQ2UID))
def split_source_map_sam(sam, out):
    cmder.run(f'split_source_map_sam {sam} {out}', msg=f'Split {sam} into small files ...', pmt=args.debug,
              stdout=sys.stdout, stderr=sys.stderr, exit_on_error=True)
    return out


@task(inputs=[], parent=split_source_map_sam, cpus=args.cpus,
      outputs=[f'{name}.mask.map.to.{STAG}.{x}{y}.target.fasta'
               for name in IDS for x in 'ATGCN' for y in 'ATGCN'])
def find_putative_target(sam, fasta):
    bam = fasta.replace('.target.fasta', '.bam')
    cmder.run(f'find_putative_target {bam} {fasta} {STAG} {args.debug}',
              msg=f'Finding putative target in {bam} ...', pmt=args.debug, debug=args.debug, exit_on_error=True)
    return fasta


@task(inputs=[], parent=find_putative_target,
      outputs=[f'{name}.mask.map.to.{STAG}.target.fasta' for name in FASTQ2UID.values()])
def merge_putative_target(path, fasta):
    name = fasta.replace('.target.fasta', '')
    cmder.run(f'cat {name}.[ATGCN][ATGCN].target.fasta > {fasta}')
    if not args.debug:
        cmder.run(f'rm {name}.[ATGCN][TGCN].target.fasta || true')
    target_csv, tmp_csv = fasta.replace(".fasta", ".csv"), fasta.replace(".fasta", ".tmp.csv")
    cmder.run(f'cat {name}.[ATGCN][ATGCN].target.csv > {tmp_csv}')
    cmder.run(f'grep -m 1 name {tmp_csv} > {target_csv}')
    cmder.run(f'grep -v name {tmp_csv} >> {target_csv}')
    cmder.run(f'rm {tmp_csv}')
    if not args.debug:
        cmder.run(f'rm {name}.[ATGCN][ATGCN].target.csv || true')
        cmder.run(f'rm {name}.[ATGCN][ATGCN].bam || true ')


@task(inputs=merge_putative_target, outputs=lambda i: i.replace('.target.fasta', '.true.target.fasta'))
def map_putative_targets_back_to_source_rna(target, true_target):
    db = bowtie2_index(f'{STAG}.bowtie2.index', args.source_rna_fasta)
    sam = target.replace('.target.fasta', '.target.back.map.sam')
    cmd = ['bowtie2', '-D', '20', '-R', '3', '-N', '0', '-L', '16', '--local', '--norc', '-i', 'S,1,0.50',
           '--score-min', 'L,12,0', '--ma', '1', '--np', '0', '--mp', '2,2', '--rdg', '5,1', '--rfg', '5,1',
           '-p', args.cpus, '--no-unal', '-x', db,
           '--un', true_target, '-f', '-U', target, '-a', '-S', sam, '2>', sam.replace('.sam', '.log')]
    cmder.run(cmd, msg=f'Mapping targets {os.path.basename(target)} back to {STAG} ...',
              pmt=args.debug, debug=args.debug, exit_on_error=True)
    return true_target


@task(inputs=[], parent=map_putative_targets_back_to_source_rna,
      outputs=[f'{name}.{STAG}.{tag}.sam'
               for name in FASTQ2UID.values() for tag in TAG])
def map_putative_target_to_target_rna(fasta, sam):
    basename, _, tag, _ = sam.rsplit('.', 3)
    db = bowtie_index(f'{tag}.bowtie.index', TARGETS[tag])
    fasta = f'{basename}.mask.map.to.{STAG}.true.target.fasta'
    cmd = ['bowtie', '-a', '--best', '--strata', '-e', 35, '-q', '-l', 8, '-n', 1, '-p', args.cpus,
           '-x', db, '--no-unal', '--norc', '-f', fasta, '--sam', sam, '--un', sam.replace('.sam', '.unmap.fasta'),
           '2>', sam.replace('.sam', '.log')]
    cmder.run(cmd, msg=f'Mapping {fasta} to {tag} ...', pmt=args.debug, debug=args.debug, exit_on_error=True)
    return sam


@task(inputs=[], parent=map_putative_target_to_target_rna, cpus=args.cpus,
      outputs=[f'{name}.{STAG}.RNA.unmap.fasta' for name in FASTQ2UID.values()])
def compile_unmapped_putative_target(inputs, fasta):
    tags = ",".join(TAG) if len(TAG) > 1 else f'{TAG[0]},'
    sample = fasta.replace(f'.{STAG}.RNA.unmap.fasta', '')
    cmder.run(f'compile_unmapped_putative_target {sample} {STAG} {tags} {fasta}',
              msg=f'Compiling unmapped putative targets for sample {sample} ...',
              pmt=args.debug, stdout=sys.stdout, stderr=sys.stderr, exit_on_error=True)


@task(inputs=compile_unmapped_putative_target, outputs=lambda i: f'{i.replace(".RNA.unmap.fasta", f".{GTAG}.sam")}')
def map_putative_target_to_genome(fasta, sam):
    prefix = sam.replace('.sam', '')
    n = sum([int(os.path.exists(f'{name}.{STAG}.{args.species}.sam')) for name in IDS])
    genome_load = 'LoadAndRemove' if n == len(IDS) - 1 else 'LoadAndKeep'
    try:
        os.makedirs(prefix, exist_ok=True)
        cmder.logger.debug(f'[ -d {prefix} ] || mkdir {prefix}')
        cmd = ['STAR', '--alignEndsType', 'EndToEnd', '--genomeDir', GENOMES[args.species],
               '--genomeLoad', genome_load, '--outBAMcompression', '10',
               '--outFileNamePrefix', f'{prefix}/',
               '--outFilterMatchNminOverLread', '0.66', '--outFilterMultimapNmax', '1',
               '--outFilterMultimapScoreRange', '1', '--outFilterScoreMin', '10',
               '--outFilterScoreMinOverLread', '0.66', '--outFilterType', 'BySJout',
               '--outReadsUnmapped', 'Fastx', '--outSAMattrRGline', 'ID:foo',
               # 'Standard' (NH HI AS nM), not 'All': this step writes SAM, and from
               # STAR 2.5 'All' pulls in the ch tag, which is BAM-only and aborts the
               # run. The published pipeline.sh for this same step used Standard too --
               # only the Python drifted. Downstream needs AS, which Standard includes.
               '--outSAMattributes', 'Standard', '--outSAMmode', 'Full',
               '--outSAMtype', 'SAM', '--outSAMunmapped', 'None',
               '--outStd', 'Log', '--readFilesIn', fasta, '--runMode', 'alignReads', '--runThreadN', args.cpus]
        cmder.run(cmd, msg=f'Mapping {fasta} to {args.species} ...', pmt=args.debug, debug=args.debug,
                  exit_on_error=True)
        cmder.run(f'{"cp" if args.debug else "mv"} {prefix}/Aligned.out.sam {sam}')
        cmder.run(f'{"cp" if args.debug else "mv"} {prefix}/Log.final.out {sam.replace(".sam", ".log")}')
    finally:
        if args.keep or args.debug:
            pass
        else:
            cmder.run(f'rm -r {prefix}')
    return sam


@task(inputs=[], cpus=args.cpus, parent=map_putative_target_to_genome,
      outputs=[f'{name}.{STAG}.{tag}.bam' for name in FASTQ2UID.values() for tag in TAGS])
def convert_sam_to_bam(sam, bam):
    sam = bam.replace('.bam', '.sam')
    if GTAG in bam and args.source_rna_bed:
        tmp_bam = bam.replace('.bam', '.tmp.bam')
        cmder.run(f'sam_to_bam {sam} {tmp_bam} True', pmt=args.debug, debug=args.debug, exit_on_error=True)
        cmder.run(f'bedtools intersect -abam {tmp_bam} -b {args.source_rna_bed} -f 0.25 -v -s > {bam}',
                  pmt=args.debug, debug=args.debug, exit_on_error=True)
        cmder.run(f'samtools index {bam}')
        if not args.debug:
            cmder.run(f'rm {tmp_bam} {tmp_bam}.bai')
    else:
        cmder.run(f'sam_to_bam {sam} {bam} True', pmt=args.debug, debug=args.debug, exit_on_error=True)
    return bam


@task(inputs=convert_sam_to_bam, outputs=lambda i: i.replace('.bam', '.chimeras.bam'))
def identify_chimeric_read(bam, out):
    cmder.run(f'identify_chimeric_read {bam} {COUNT} {out}', exit_on_error=True,
              msg=f'Identifying chimeric reads in {bam}', pmt=args.debug, debug=args.debug)


@task(inputs=identify_chimeric_read, outputs=lambda i: i.replace('.chimeras.bam', '.chimeras.pos.bw'))
def chimeric_bigwig(bam, bw):
    cmder.run(f'bam_to_bw {bam}', msg=f'Creating BigWig tracks using {bam} ...', pmt=args.debug,
              debug=args.debug, exit_on_error=True)


# @task(inputs=[f'{uid}.{STAG}.{args.species}.chimeras.bam' for uid in IDS], parent=identify_chimeric_read,
#       outputs=lambda i: i.replace('.chimeras.bam', '.chimeras.peak.bed'))
# def genomic_chimeric_peak_cluster(bam, bed):
#     if not os.path.exists(bed):
#         cmd = f'clipper --species {args.species} --processors {args.cpus} --bam {bam} --outfile {bed}'
#         cmder.run(cmd, msg=f'Calling peaks in {bam} ...', debug=args.debug, pmt=args.debug, exit_on_error=True)


@task(inputs=[], parent=identify_chimeric_read, outputs=['individual.chimeras.count.tsv'])
def individual_chimeras(inputs, outputs):
    for uid in IDS:
        for tag in TAGS:
            bam = f'{uid}.{STAG}.{tag}.chimeras.bam'
            cmd = ['parse_individual_chimeras', bam,
                   '--uid', uid, '--tag', tag, '--cpus', args.cpus,
                   '--minimum_chimeras', args.minimum_chimeras, '--count_output', outputs]
            if tag == GTAG:
                if args.clipper_rna_species:
                    cmd.extend(['--species', GTAG])
            elif tag == args.target_rRNA_tag and args.clipper_rna_species:
                cmd.extend(['--species', args.clipper_rna_species])
            cmder.run(cmd, msg=f'Parsing {tag} individual chimeras for {uid} ...', debug=True, exit_on_error=True)


def cleanup():
    log = ' '.join([f'{uid}*.log' for uid in IDS if os.path.exists(f'{uid}*.log')])
    if log:
        cmder.run(f'tail -n +1 {log} > map.metrics.log')
        cmder.run(f'rm {log}')
    
    cmder.run(f"rm {' '.join([f'{uid}.*.trim.fastq.gz' for uid in IDS])} || true")
    cmder.run(f"rm {' '.join([f'{uid}.*.target.*' for uid in IDS])} || true")
    cmder.run(f"rm {' '.join([f'{uid}.*.fastq' for uid in IDS])} || true")
    cmder.run(f"rm {' '.join([f'{uid}.*.mask.fasta' for uid in IDS])} || true")
    cmder.run(f"rm {' '.join([f'{uid}.*.duplicated.csv' for uid in IDS])} || true")
    cmder.run(f"rm {' '.join([f'{uid}.*.sam' for uid in IDS])} || true")
    cmder.run(f"rm {' '.join([f'{uid}.[ATGCN][ATGCN].bam' for uid in IDS])} || true")
    cmder.run(f"rm {' '.join([f'{uid}.candidate.bam*' for uid in IDS])} || true")
    cmder.run(f"rm {' '.join([f'{uid}.unmap.fasta' for uid in IDS])} || true")
    cmder.run(f"rm -r {' '.join([f'{uid}*.repeat.map' for uid in IDS])} || true")
    cmder.run(f"rm -r {' '.join([f'{uid}*.genome.map' for uid in IDS])} || true")
    cmder.run(f"rm -r {' '.join([f'{uid}*.genome.length.tsv' for uid in IDS])} || true")
    cmder.run(f"rm -r {' '.join([f'{uid}*.bw.*.bg' for uid in IDS])} || true")
    cmder.run(f"rm -r {' '.join([f'{uid}*.{STAG}.{GTAG}' for uid in IDS])} || true")
    cmder.run(f"tail -n +1 {' '.join([f'{uid}*.log' for uid in IDS])} > map.metric.and.log")
    cmder.run(f"rm {' '.join([f'{uid}*.log' for uid in IDS])} || true")
    cmder.run(f"rm -f {STAG}*bt2 || true")
    cmder.run(f"rm -r {' '.join([f'{tag}*.ebwt' for tag in TAG])} || true")
    

@task(inputs=[], outputs=[f'{".".join(IDS)}.{STAG}.{".".join(TAGS)}.chimeras.report.html'],
      parent=identify_chimeric_read)
def report(inputs, outputs):
    tags, uid = ','.join(TAGS), ','.join(IDS) if len(IDS) > 1 else f'{IDS[0]},'
    cmder.run(f'sno_chimeras_summary --stag {STAG} --tags {tags} --ids {uid} '
              f'--output {outputs} --count {COUNT} --gtag {GTAG}',
              msg='Generating summary report ...', exit_on_error=True,
              pmt=args.debug, stdout=sys.stdout, stderr=sys.stderr)
    if args.keep or args.debug:
        pass
    else:
        logger.info('Cleaning up ...')
        cleanup()
    logger.info('Mission accomplished.')
    
  
@logger.catch()
def main():
    flow = Flow(name='miR Chimeras', description='Identify miR chimeras using bowtie and STAR.')
    try:
        flow.run(cpus=args.cpus, dry_run=args.dryrun)
    finally:
        cmd = """for i in $(ipcs -m | grep "$(whoami)" | awk '{ print $2 }'); do ipcrm -m "$i"; done"""
        cmder.run(cmd, executable='/usr/bin/bash', fmt_cmd=False, log_cmd=False)
    

if __name__ == '__main__':
    main()
