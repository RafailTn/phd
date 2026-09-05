#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Compare one pipeline arm against the published hg19 chimeras, read name by read name.

Read names carry the original Illumina identifier plus the UMI and are therefore
genome-independent, so an hg38 arm can be compared to the hg19 publication
directly -- what differs is which reads become chimeras, not what they are called.

For every published chimera the arm does not reproduce, the read name is traced
through the arm's own kept intermediates to say *where* it was lost:

    <uid>.mask.fasta                                   survived repeat+genome masking
    <uid>.mask.map.to.<stag>.target.fasta              guide hit, target arm extracted
    <uid>.mask.map.to.<stag>.true.target.fasta         survived the back-map filter

anything still present after the last of those was lost at the genome-mapping
step itself, which is the only stage that depends on the genome build.

The target-arm FASTA headers carry a trailing '_<offset>' the pipeline appends
when it cuts the arm out of the read, so those names are truncated at the last
underscore before comparison; the mask FASTA and the CSVs use the bare name.
"""

import argparse
import csv
import os
import sys
from statistics import median


def published_reads(path):
    """{read_name: target arm length} from the published chimeras CSV."""
    reads = {}
    with open(path, newline='') as f:
        for row in csv.DictReader(f):
            try:
                reads[row['read_name']] = int(row['map_to_hg19_length'])
            except (KeyError, ValueError):
                reads[row['read_name']] = 0
    return reads


def arm_reads(path):
    with open(path, newline='') as f:
        return {row['read_name'] for row in csv.DictReader(f)}


def fasta_names(path, strip_offset=False):
    names = set()
    if not os.path.exists(path):
        return None
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()[0].split('/')[0]
                if strip_offset:
                    name = name.rsplit('_', 1)[0]
                names.add(name)
    return names


def pct(n, d):
    return f'{100 * n / d:.1f}%' if d else 'n/a'


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--outdir', required=True, help='Pipeline output directory for the arm.')
    p.add_argument('--uid', required=True, help='Sample identifier.')
    p.add_argument('--gtag', default='hg38', help='Genome tag of the arm, default: %(default)s.')
    p.add_argument('--stag', default='snoRNA', help='Source RNA tag, default: %(default)s.')
    p.add_argument('--published', required=True, help='Published hg19 chimeras CSV.')
    p.add_argument('--other-tags', default='rRNA,snRNA,tRNA',
                   help='Target classes tried before the genome, comma separated. '
                        'A read the pipeline assigned to one of these never reached '
                        'the genome step, so it is reported separately rather than '
                        'counted as an alignment failure. Default: %(default)s.')
    p.add_argument('--label', default='', help='Name for this arm in the output.')
    p.add_argument('--out', default='', help='Write the summary here as well as to stdout.')
    a = p.parse_args()

    arm_csv = os.path.join(a.outdir, f'{a.uid}.{a.stag}.{a.gtag}.chimeras.csv')
    if not os.path.exists(arm_csv):
        sys.exit(f'no chimeras CSV at {arm_csv}')

    pub = published_reads(a.published)
    ours = arm_reads(arm_csv)
    shared, missed, extra = set(pub) & ours, set(pub) - ours, ours - set(pub)

    lines = [f'# {a.label or a.outdir} vs published hg19', '',
             f'{"published hg19 chimeras":34s} {len(pub):>8,}',
             f'{"this arm":34s} {len(ours):>8,}',
             f'{"shared read names":34s} {len(shared):>8,}'
             f'   ({pct(len(shared), len(pub))} of published, {pct(len(shared), len(ours))} of this arm)',
             f'{"published only":34s} {len(missed):>8,}',
             f'{"this arm only":34s} {len(extra):>8,}', '']

    stages = [
        ('survived repeat+genome masking',
         fasta_names(os.path.join(a.outdir, f'{a.uid}.mask.fasta'))),
        ('guide hit, target arm extracted',
         fasta_names(os.path.join(a.outdir, f'{a.uid}.mask.map.to.{a.stag}.target.fasta'), True)),
        ('survived back-map filter',
         fasta_names(os.path.join(a.outdir, f'{a.uid}.mask.map.to.{a.stag}.true.target.fasta'), True)),
    ]
    if all(s is not None for _, s in stages):
        lines.append(f'where the {len(missed):,} published-only reads were lost:')
        remaining = missed
        for label, names in stages:
            kept = remaining & names
            lines.append(f'  {label:36s} kept {len(kept):>7,}   lost {len(remaining) - len(kept):>7,}')
            remaining = kept
        # "Reached the genome step" is where the residual used to stop, but the
        # pipeline tries rRNA/snRNA/tRNA *before* the genome and first hit wins.
        # A read captured there was still called a chimera -- just with a
        # different target -- so counting it as an alignment failure overstates
        # the disagreement. Split the residual out.
        elsewhere = set()
        for tag in [t for t in a.other_tags.split(',') if t]:
            other = os.path.join(a.outdir, f'{a.uid}.{a.stag}.{tag}.chimeras.csv')
            if os.path.exists(other):
                elsewhere |= (remaining & arm_reads(other))
        unplaced = remaining - elsewhere
        lines.append(f'  {"-> reached the genome step":36s}      {len(remaining):>7,}'
                     f'   ({pct(len(remaining), len(missed))} of the loss)')
        if elsewhere:
            lines.append('')
            lines.append(f'of those {len(remaining):,}, this arm did call a chimera, '
                         f'against a different target:')
            for tag in [t for t in a.other_tags.split(',') if t]:
                other = os.path.join(a.outdir, f'{a.uid}.{a.stag}.{tag}.chimeras.csv')
                if os.path.exists(other):
                    n = len(remaining & arm_reads(other))
                    if n:
                        lines.append(f'  {"claimed by " + tag:36s}      {n:>7,}'
                                     f'   ({pct(n, len(missed))} of the loss)')
            lines.append(f'  {"genuinely unplaced":36s}      {len(unplaced):>7,}'
                         f'   ({pct(len(unplaced), len(missed))} of the loss)')
    else:
        lines.append('intermediates not kept (--keep), so the loss cannot be attributed to a stage')
    lines.append('')

    for label, group in (('published-only', missed), ('shared', shared)):
        lengths = sorted(pub[n] for n in group)
        if lengths:
            short = sum(1 for x in lengths if x <= 25)
            lines.append(f'{label + " target-arm length":34s} n={len(lengths):>7,}  '
                         f'median {median(lengths):>4.0f} nt   <=25 nt {pct(short, len(lengths))}')

    text = '\n'.join(lines)
    print(text)
    if a.out:
        with open(a.out, 'w') as o:
            o.write(text + '\n')


if __name__ == '__main__':
    main()
