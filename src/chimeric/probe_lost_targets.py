#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""What does --outFilterMultimapNmax 1 throw away, and does it fall on AluACAs?

The genome step keeps a target arm only if it places at exactly one locus. An arm
inside an Alu almost never does, and an antisense Alu is the partner an Alu-derived
guide is most likely to pair with -- so the filter may be removing the very class the
AluACA hypothesis is about, before anything is counted. The reports measure AluACA
pairing on whatever survives; this measures what does not.

Every putative target arm entering the genome step is re-mapped with the same settings
except --outFilterMultimapNmax, opened from 1 to `--nmax`. Each arm then falls into:

  unique      NH == 1, the pipeline keeps it
  multi       NH > 1, the pipeline discards it with no record -- what this is about
  unmapped    never reached the multimapping decision (the 0.66 thresholds)

Multi arms are split by what their loci are: an Alu in the same orientation as the arm
(cannot base-pair with an Alu-derived guide), an antisense Alu (can), or no Alu at all.
Guide class comes from the pipeline's own putative-target table, so arms that never
became calls still carry one, and input and IP are measured identically.

    python3 src/chimeric/probe_lost_targets.py \\
        --outdir results/chimeric_chr25/arm3_hg19_merged/SRR30692552 --uid SRR30692552 \\
        --index ref/chimeric/hg19_chr25_star_index --rmsk ref/chimeric/rmsk.hg19.bed \\
        --out results/chimeric_chr25/SRR30692552.lost_targets.tsv
"""

import argparse
import csv
import os
import subprocess
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import find_input, find_tool
from annotate_chimeras import load_alu_names, guide_class
from probe_multimapping import run_star, base_name, ref_span


def read_target_table(path, stag):
    """{read_name: guide names} from the pipeline's putative-target CSV.

    Column names come from upstream and have changed between versions, so the read-name
    and guide columns are found by shape rather than assumed.
    """
    with open(path) as fh:
        head = next(csv.reader(fh))
    name_col = next((c for c in head if c.lower() in ('name', 'read_name', 'qname')), None)
    guide_col = next((c for c in head if c == f'reference_{stag}'), None)
    if not (name_col and guide_col):
        sys.exit(f'{path}: need a read-name column and reference_{stag}; found {head}')
    # Read names end with _<UMI>, and the FASTA headers add _<offset> on top, so the two
    # files can key on either form. Store both and let the lookup decide.
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            name = row[name_col]
            out[name] = row[guide_col]
            tail = name.rsplit('_', 1)
            if len(tail) == 2 and tail[1].isdigit():
                out[tail[0]] = row[guide_col]
    return out


def alu_loci(sam_hits, rmsk, bedtools, workdir):
    """{read: set of 'sense'/'antisense'/'none'} for the loci of each multimapping arm."""
    os.makedirs(workdir, exist_ok=True)
    bed = os.path.join(workdir, 'loci.bed')
    keys = []
    with open(bed, 'w') as o:
        for read, hits in sam_hits.items():
            for i, (chrom, pos, span, nh, score, strand) in enumerate(hits):
                if nh <= 1:
                    continue
                keys.append(read)
                o.write(f'{chrom}\t{pos - 1}\t{pos - 1 + span}\t{len(keys) - 1}\t0\t{strand}\n')
    if not keys:
        return {}
    subprocess.run(f'LC_ALL=C sort -k1,1 -k2,2n {bed} -o {bed}', shell=True, check=True)
    alu = os.path.join(workdir, 'alu.bed')
    if not os.path.exists(alu):
        subprocess.run(f"awk -F'\\t' '$4 ~ /^Alu/' {rmsk} | LC_ALL=C sort -k1,1 -k2,2n > {alu}",
                       shell=True, check=True)
    out = defaultdict(set)
    for orient, flag in (('sense', '-s'), ('antisense', '-S')):
        q = subprocess.run(f'{bedtools} intersect -a {bed} -b {alu} -u {flag} -f 0.5 | cut -f4',
                           shell=True, capture_output=True, text=True)
        for line in q.stdout.split():
            out[keys[int(line)]].add(orient)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--outdir', required=True, help='Pipeline output directory for the sample.')
    p.add_argument('--uid', required=True)
    p.add_argument('--stag', default='snoRNA')
    p.add_argument('--index', default=os.environ.get('GENOME_INDEX', ''),
                   help='STAR index the pipeline places target arms on.')
    p.add_argument('--rmsk', default=os.environ.get('RMSK_BED', ''), help='RepeatMasker BED.')
    p.add_argument('--alu-fasta', default=None, help='AluACA union FASTA, to class guides.')
    p.add_argument('--nmax', type=int, default=100,
                   help='--outFilterMultimapNmax for the probe, default: %(default)s.')
    p.add_argument('--work', default='', help='Scratch directory for the STAR run.')
    p.add_argument('--star', default=None)
    p.add_argument('--bedtools', default=None)
    p.add_argument('--cpus', type=int, default=int(os.environ.get('CPUS') or os.cpu_count() or 1))
    p.add_argument('--out', required=True, help='Per-read TSV to write.')
    a = p.parse_args()
    a.star = find_tool('STAR', a.star)
    a.bedtools = find_tool('bedtools', a.bedtools)
    a.alu_fasta = a.alu_fasta or find_input('AluACA_union_nr.fasta')
    work = a.work or os.path.join(a.outdir, 'lost_targets_work')

    fasta = os.path.join(a.outdir, f'{a.uid}.{a.stag}.RNA.unmap.fasta')
    table = os.path.join(a.outdir, f'{a.uid}.mask.map.to.{a.stag}.target.csv')
    for f in (fasta, table):
        if not os.path.exists(f):
            sys.exit(f'missing {f}')
    alu_names = load_alu_names(a.alu_fasta)
    guides = read_target_table(table, a.stag)
    arms = [l[1:].strip() for l in open(fasta) if l.startswith('>')]
    reads = {base_name(h) for h in arms}
    print(f'{len(arms):,} target arms from {len(reads):,} reads entering the genome step')

    sam = run_star(a.star, a.index, fasta, os.path.join(work, 'star'), a.cpus, a.nmax, 1)
    hits = defaultdict(list)
    with open(sam) as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            f = line.rstrip('\n').split('\t')
            tags = {t.split(':', 2)[0]: t.split(':', 2)[2] for t in f[11:]}
            hits[base_name(f[0])].append((f[2], int(f[3]), ref_span(f[5]),
                                          int(tags.get('NH', 1)), int(tags.get('AS', 0)),
                                          '-' if int(f[1]) & 16 else '+'))
    print(f'  {len(hits):,} reads placed somewhere at --outFilterMultimapNmax {a.nmax}')
    orient = alu_loci(hits, a.rmsk, a.bedtools, work)

    rows, tally = [], defaultdict(Counter)
    for read in reads:
        g = guides.get(read)
        cls = guide_class(g, alu_names) if g else 'unknown'
        h = hits.get(read)
        if not h:
            fate, loci = 'unmapped', 0
        else:
            loci = max(x[3] for x in h)
            fate = 'unique (kept)' if loci == 1 else 'multi (discarded)'
        o = orient.get(read, set())
        alu = ('antisense Alu' if 'antisense' in o else
               'sense Alu only' if 'sense' in o else 'no Alu locus')
        rows.append((read, cls, fate, loci, alu if fate.startswith('multi') else ''))
        tally[cls][fate] += 1
        if fate.startswith('multi'):
            tally[cls][f'  of those: {alu}'] += 1
    with open(a.out, 'w') as o:
        o.write('read_name\tguide_class\tfate\tloci\tmulti_alu\n')
        o.writelines('\t'.join(map(str, r)) + '\n' for r in rows)
    print(f'\nwrote {a.out}')

    print('\n=== what the genome step keeps and discards, by guide class ===')
    for cls in sorted(tally):
        n = sum(v for k, v in tally[cls].items() if not k.startswith('  '))
        print(f'\n{cls}  ({n:,} reads)')
        for k in ('unique (kept)', 'multi (discarded)', 'unmapped'):
            print(f'  {k:<22} {tally[cls][k]:>8,}  {100.0 * tally[cls][k] / max(n, 1):5.1f}%')
        m = tally[cls]['multi (discarded)']
        for k in ('  of those: antisense Alu', '  of those: sense Alu only', '  of those: no Alu locus'):
            print(f'  {k:<22} {tally[cls][k]:>8,}  {100.0 * tally[cls][k] / max(m, 1):5.1f}% of discarded')


if __name__ == '__main__':
    main()
