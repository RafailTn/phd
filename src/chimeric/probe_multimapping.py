#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Why does STAR discard the published target arms this arm fails to place?

The genome step runs --outFilterMultimapNmax 1, so a read that matches two loci
is thrown away with no record of where it matched. Log.final.out only gives
bucket totals: for arm1 that is 39,855 reads "mapped to too many loci" and
52,102 "unmapped: too short", and the 12,505 published calls we miss could sit
in either. Totals cannot say which.

This re-maps those arms against the same index with the filter opened to 20 and
NH reported, so each read states its own reason. Three outcomes, three culprits:

  NH 2-3 with a hit on an unplaced contig   the index carries scaffolds the
                                            published one (25 references) does
                                            not -> rebuild with 25 contigs;
  NH large, every hit a main chromosome     no reference change helps; this is
                                            how 2.7.11b enumerates loci where
                                            2.4.0j stopped -> STAR version;
  still unmapped at NH<=20                  never reached the multimapping
                                            decision -> the 0.66 length/score
                                            thresholds, not the index.

Reads we *do* share with the publication are profiled alongside as a control:
without them a scaffold hit among the missed reads means nothing, because the
shared reads may carry them at the same rate.

    python3 src/chimeric/probe_multimapping.py \\
        --outdir results/chimeric/arm1_hg19_plain/SRR30692552 \\
        --uid SRR30692552 --gtag hg19 \\
        --index ref/chimeric/hg19_star_index \\
        --published data/DKC1_IP.snoRNA.hg19.chimeras.csv
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import find_input, find_tool, require

# The published output uses exactly these 25 references and no scaffolds, which
# is what makes "is this hit on a main chromosome?" the discriminating question.
MAIN = re.compile(r'^chr(\d{1,2}|X|Y|M)$')


def read_fasta(path):
    """[(header, sequence)] in file order. Headers are '<read_name>_<offset>',
    so several records can share a read name; all of them are kept."""
    out, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    out.append((name, ''.join(buf)))
                name, buf = line[1:].split()[0].split('/')[0], []
            else:
                buf.append(line.strip())
    if name:
        out.append((name, ''.join(buf)))
    return out


def base_name(header):
    """'<read_name>_<offset>' -> '<read_name>'. The offset is appended by
    find_putative_target and is not part of the name the publication uses."""
    return header.rsplit('_', 1)[0]


def published_loci(path):
    """{read_name: (chrom, start, stop)} for the published genomic placement."""
    loci = {}
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh):
            try:
                loci[row['read_name']] = (row['reference_hg19'],
                                          int(row['map_to_hg19_ref_start']),
                                          int(row['map_to_hg19_ref_stop']))
            except (KeyError, ValueError):
                loci[row['read_name']] = None
    return loci


def arm_calls(path):
    with open(path, newline='') as fh:
        return {row['read_name'] for row in csv.DictReader(fh)}


def write_subset(records, wanted, path):
    """Write every record whose read name is in `wanted`; return the names hit."""
    seen = set()
    with open(path, 'w') as fh:
        for header, seq in records:
            if base_name(header) in wanted:
                fh.write(f'>{header}\n{seq}\n')
                seen.add(base_name(header))
    return seen


def run_star(star, index, fasta, prefix, cpus, nmax, score_range):
    """The pipeline's genome step (sno-chimeras.py:405-418) with two changes:
    --outFilterMultimapNmax opened from 1 to `nmax`, and NH/HI added to the SAM
    attributes so the count of loci survives into the output. Everything else --
    EndToEnd, the 0.66 thresholds, BySJout -- is held fixed, or the run would
    not be measuring the same decision the pipeline makes.

    `score_range` is --outFilterMultimapScoreRange, 1 in the pipeline. At 1 an
    alignment two points worse than the best is not reported at all, so a read
    matching the rDNA scaffold GL000220.1 can hide a real chromosomal alignment
    that would become the best one once the scaffold is gone. Widening it says
    whether the alignment exists."""
    os.makedirs(prefix, exist_ok=True)
    cmd = [star, '--alignEndsType', 'EndToEnd', '--genomeDir', index,
           '--genomeLoad', 'NoSharedMemory',
           '--outFileNamePrefix', f'{prefix}/',
           '--outFilterMatchNminOverLread', '0.66',
           '--outFilterMultimapNmax', str(nmax),
           '--outFilterMultimapScoreRange', str(score_range),
           '--outFilterScoreMin', '10',
           '--outFilterScoreMinOverLread', '0.66',
           '--outFilterType', 'BySJout',
           '--outReadsUnmapped', 'Fastx',
           '--outSAMattributes', 'NH', 'HI', 'AS', 'nM',
           '--outSAMmode', 'Full', '--outSAMtype', 'SAM',
           '--outSAMunmapped', 'None', '--outStd', 'Log',
           '--readFilesIn', fasta, '--runMode', 'alignReads',
           '--runThreadN', str(cpus)]
    print('  ' + ' '.join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)
    return os.path.join(prefix, 'Aligned.out.sam')


_CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')


def ref_span(cigar):
    """Bases of the reference the alignment covers."""
    return sum(int(n) for n, op in _CIGAR.findall(cigar) if op in 'MDN=X')


def parse_sam(path):
    """{read_name: [(chrom, pos, span, nh)]}, pooled over the several arms a
    read can contribute (headers carry the offset, read names do not)."""
    hits = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            f = line.rstrip('\n').split('\t')
            nh = 1
            for tag in f[11:]:
                if tag.startswith('NH:i:'):
                    nh = int(tag[5:])
                    break
            score = 0
            for tag in f[11:]:
                if tag.startswith('AS:i:'):
                    score = int(tag[5:])
                    break
            hits.setdefault(base_name(f[0]), []).append(
                (f[2], int(f[3]), ref_span(f[5]), nh, score))
    return hits


def overlaps(chrom, pos, span, locus, slack=10):
    if not locus:
        return False
    c, s, e = locus
    if c != chrom:
        return False
    lo, hi = min(s, e) - slack, max(s, e) + slack
    return pos <= hi and pos + span >= lo


def simulate_main_only(label, names, hits, loci, check_locus, score_range=1):
    """What the 25-reference index would have done, computed from this SAM.

    Discard every alignment off the main chromosomes, then apply the pipeline's
    own rule to what is left: keep alignments within `score_range` of the best,
    and require exactly one. This is a prediction, not a rerun -- removing
    sequence also changes STAR's seed search and its BySJout pass -- but it uses
    the real alignments and scores, so it bounds what the rebuild can deliver.
    """
    unique = at_published = still_multi = gone = 0
    for n in names:
        main = [a for a in hits.get(n, []) if MAIN.match(a[0])]
        if not main:
            gone += 1
            continue
        best = max(a[4] for a in main)
        keep = [a for a in main if a[4] >= best - score_range]
        if len(keep) > 1:
            still_multi += 1
            continue
        unique += 1
        if check_locus and overlaps(keep[0][0], keep[0][1], keep[0][2], loci.get(n)):
            at_published += 1

    total = len(names) or 1
    print(f'\n{label}: predicted outcome with the scaffolds removed')
    for lbl, v in (('uniquely mapped', unique),
                   ('still multimapping', still_multi),
                   ('no alignment left', gone)):
        print(f'  {lbl:<22} {v:>8,}  {100.0 * v / total:5.1f}%')
    if check_locus:
        print(f'  {"of those, at the published locus":<22} '
              f'{at_published:>8,}  {100.0 * at_published / total:5.1f}%')


def profile(label, names, hits, loci, check_locus):
    """One block of the report: where the reads went once the filter was opened."""
    nh_bucket = Counter()
    unmapped = 0
    with_scaffold = 0
    with_main = 0
    at_published = 0
    for n in sorted(names):
        al = hits.get(n)
        if not al:
            unmapped += 1
            nh_bucket['unmapped'] += 1
            continue
        nh = max(a[3] for a in al)
        nh_bucket['1' if nh == 1 else '2-3' if nh <= 3 else '4-10' if nh <= 10
                  else '11-20'] += 1
        if any(not MAIN.match(a[0]) for a in al):
            with_scaffold += 1
        if any(MAIN.match(a[0]) for a in al):
            with_main += 1
        if check_locus and any(overlaps(a[0], a[1], a[2], loci.get(n)) for a in al):
            at_published += 1

    total = len(names) or 1
    print(f'\n{label}  n = {len(names):,}')
    for k in ('1', '2-3', '4-10', '11-20', 'unmapped'):
        v = nh_bucket[k]
        tag = 'loci' if k != 'unmapped' else ''
        print(f'  {k:>9} {tag:<5} {v:>8,}  {100.0 * v / total:5.1f}%')
    mapped = len(names) - unmapped
    if mapped:
        print(f'  of the {mapped:,} that mapped at all:')
        print(f'    at least one hit off the 25 main chromosomes '
              f'{with_scaffold:>7,}  {100.0 * with_scaffold / mapped:5.1f}%')
        print(f'    at least one hit ON a main chromosome        '
              f'{with_main:>7,}  {100.0 * with_main / mapped:5.1f}%')
        if check_locus:
            print(f'    published locus among the hits              '
                  f'{at_published:>7,}  {100.0 * at_published / mapped:5.1f}%')


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--outdir', required=True, help='Pipeline output directory for the arm.')
    p.add_argument('--uid', required=True, help='Sample identifier.')
    p.add_argument('--gtag', default='hg19', help='Genome tag of the arm, default: %(default)s.')
    p.add_argument('--stag', default='snoRNA', help='Source RNA tag, default: %(default)s.')
    p.add_argument('--index', default='', help='STAR genome index; default $GENOME_INDEX.')
    p.add_argument('--published', default='', help='Published hg19 chimeras CSV.')
    p.add_argument('--work', default='', help='Where to put the probe fastas and STAR output.')
    p.add_argument('--nmax', type=int, default=20, help='--outFilterMultimapNmax, default: %(default)s.')
    p.add_argument('--score-range', type=int, default=1,
                   help='--outFilterMultimapScoreRange, default: %(default)s (the '
                        'pipeline value). Widen it to reveal alignments the rDNA '
                        'scaffold outscores.')
    p.add_argument('--cpus', type=int, default=int(os.environ.get('CPUS') or os.cpu_count() or 4))
    p.add_argument('--star', default='', help='STAR executable.')
    p.add_argument('--from-sam', default='',
                   help='Skip STAR and read alignments from an existing run: pass '
                        'a path containing the literal MISSED, which is replaced by '
                        'missed/shared, e.g. work/multimap_probe_wide/SRR30692552/'
                        'MISSED/Aligned.out.sam.')
    a = p.parse_args()

    published = a.published or find_input('DKC1_IP.snoRNA.hg19.chimeras.csv')
    index = a.index or os.environ.get('GENOME_INDEX') or ''
    star = find_tool('STAR', a.star)
    unmap = os.path.join(a.outdir, f'{a.uid}.{a.stag}.RNA.unmap.fasta')
    calls = os.path.join(a.outdir, f'{a.uid}.{a.stag}.{a.gtag}.chimeras.csv')
    need = [('published chimeras', published), ('genome-step input', unmap),
            ('arm chimeras', calls)]
    if not a.from_sam:      # re-reading a finished run needs neither
        need += [('STAR index', index), ('STAR', star)]
    require(need)

    work = a.work or os.path.join(os.environ.get('WORK') or 'work',
                                  'multimap_probe', os.path.basename(a.outdir))
    os.makedirs(work, exist_ok=True)

    loci = published_loci(published)
    ours = arm_calls(calls)
    records = read_fasta(unmap)
    reached = {base_name(h) for h, _ in records}

    # Only reads that actually reached the genome step can have been discarded
    # by it; anything lost earlier is a different question, already answered by
    # compare_to_published.py.
    missed = (set(loci) - ours) & reached
    shared = (set(loci) & ours) & reached

    print(f'published chimeras            {len(loci):>8,}')
    print(f'called by this arm            {len(ours):>8,}')
    print(f'published-only, reached STAR  {len(missed):>8,}')
    print(f'shared, reached STAR          {len(shared):>8,}   (control)')

    check_locus = (a.gtag == 'hg19')
    if not check_locus:
        print(f'\n!!! --gtag {a.gtag}: published coordinates are hg19, so the '
              f'"published locus among the hits" check is skipped.')

    for label, names in (('MISSED (published-only)', missed),
                         ('SHARED (control)', shared)):
        tag = 'missed' if 'MISSED' in label else 'shared'
        fasta = os.path.join(work, f'{tag}.fasta')
        hit_names = write_subset(records, names, fasta)
        sam = (a.from_sam.replace('MISSED', tag) if a.from_sam
               else run_star(star, index, fasta, os.path.join(work, tag),
                             a.cpus, a.nmax, a.score_range))
        hits = parse_sam(sam)
        profile(label, hit_names, hits, loci, check_locus)
        simulate_main_only(label, hit_names, hits, loci, check_locus)

    print(f'\nSTAR output kept under {work}')


if __name__ == '__main__':
    main()
