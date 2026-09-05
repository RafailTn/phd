#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Why did published chimeras this arm reached the genome step with fail to place?

`compare_to_published.py` narrows the disagreement to reads that survived masking,
got a guide hit, passed the back-map filter, were not claimed by another target
class, and still produced no genomic alignment. That is a residual, not a cause.

There are only two ways to get there:

  1. the arm we handed STAR is not the arm the publication placed -- the guide
     match landed differently, so the read was split at a different point, or the
     other flank was chosen (find_putative_target keeps only the longer flank);
  2. the arm is the same and STAR did not place it -- an index or parameter
     difference (EndToEnd, --outFilterMultimapNmax 1, the 0.66 thresholds).

Those need opposite fixes, so this script separates them by reconstructing the
published target arm from the published CSV and comparing it to ours, base for base.

    python3 src/chimeric/diagnose_unplaced.py \\
        --outdir results/chimeric/arm1_hg19_plain/SRR30692552 \\
        --uid SRR30692552 --gtag hg19 \\
        --published data/DKC1_IP.snoRNA.hg19.chimeras.csv
"""

import argparse
import csv
import os
import sys
from statistics import median


def read_fasta(path):
    """{name: sequence}. Headers are '<read_name>_<offset>'; the offset is where
    the arm starts in the read, so it is kept separately rather than discarded."""
    out = {}
    if not os.path.exists(path):
        return out
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    out[name] = ''.join(buf)
                name, buf = line[1:].split()[0].split('/')[0], []
            else:
                buf.append(line.strip())
    if name:
        out[name] = ''.join(buf)
    return out


def split_offset(header):
    base, _, off = header.rpartition('_')
    try:
        return base, int(off)
    except ValueError:
        return header, None


def csv_names(path):
    if not os.path.exists(path):
        return set()
    with open(path, newline='') as fh:
        return {r['read_name'] for r in csv.DictReader(fh)}


def pct(n, d):
    return f'{100 * n / d:.1f}%' if d else 'n/a'


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--outdir', required=True, help='Pipeline output directory for the arm.')
    p.add_argument('--uid', required=True, help='Sample identifier.')
    p.add_argument('--stag', default='snoRNA', help='Source RNA tag, default: %(default)s.')
    p.add_argument('--gtag', default='hg38', help='Genome tag, default: %(default)s.')
    p.add_argument('--published', required=True, help='Published chimeras CSV.')
    p.add_argument('--examples', type=int, default=5, help='Worked examples to print.')
    p.add_argument('--out', default='', help='Write the summary here as well as to stdout.')
    a = p.parse_args()

    pub = {}
    with open(a.published, newline='') as fh:
        for r in csv.DictReader(fh):
            try:
                s, e = int(r['map_to_hg19_read_start']), int(r['map_to_hg19_read_stop'])
                gs, ge = int(r['map_to_snoRNA_read_start']), int(r['map_to_snoRNA_read_stop'])
            except (KeyError, ValueError):
                continue
            pub[r['read_name']] = dict(seq=r['sequence'], t0=s, t1=e, g0=gs, g1=ge)

    ours = csv_names(os.path.join(a.outdir, f'{a.uid}.{a.stag}.{a.gtag}.chimeras.csv'))
    missed = set(pub) - ours

    arms_path = os.path.join(a.outdir, f'{a.uid}.mask.map.to.{a.stag}.true.target.fasta')
    raw = read_fasta(arms_path)
    if not raw:
        sys.exit(f'no target arms at {arms_path} (the run needs --keep)')
    arms = {}
    for header, seq in raw.items():
        base, off = split_offset(header)
        arms[base] = (seq, off)

    # Only reads we actually handed to the genome step can be diagnosed here.
    subject = [n for n in missed if n in arms]

    same, diff, shorter, longer = [], [], [], []
    for n in subject:
        our_seq = arms[n][0]
        # The published target arm, reconstructed from its own read coordinates.
        their_seq = pub[n]['seq'][pub[n]['t0']:pub[n]['t1']]
        if our_seq == their_seq:
            same.append(n)
        else:
            diff.append(n)
            if len(our_seq) < len(their_seq):
                shorter.append(n)
            elif len(our_seq) > len(their_seq):
                longer.append(n)

    L = []
    L.append(f'# why {len(subject):,} unplaced published chimeras failed  ({a.outdir})')
    L.append('')
    L.append(f'{"published chimeras":38s} {len(pub):>8,}')
    L.append(f'{"not called by this arm":38s} {len(missed):>8,}')
    L.append(f'{"...with a target arm we can inspect":38s} {len(subject):>8,}')
    L.append('')
    L.append('did we hand STAR the same sequence the publication placed?')
    L.append(f'  {"identical arm -> STAR/index difference":38s} {len(same):>8,}   '
             f'({pct(len(same), len(subject))})')
    L.append(f'  {"different arm -> split-point difference":38s} {len(diff):>8,}   '
             f'({pct(len(diff), len(subject))})')
    if diff:
        L.append(f'     {"ours shorter than theirs":35s} {len(shorter):>8,}')
        L.append(f'     {"ours longer than theirs":35s} {len(longer):>8,}')
        ol = [len(arms[n][0]) for n in diff]
        tl = [len(pub[n]['seq'][pub[n]['t0']:pub[n]['t1']]) for n in diff]
        L.append(f'     {"our arm length, median":35s} {median(ol):>8.0f} nt')
        L.append(f'     {"their arm length, median":35s} {median(tl):>8.0f} nt')
        # An arm this short cannot survive --outFilterMultimapNmax 1 in a 3 Gb genome.
        tiny = sum(1 for x in ol if x < 25)
        L.append(f'     {"our arm < 25 nt (unmappable)":35s} {tiny:>8,}   ({pct(tiny, len(diff))})')
    if same:
        sl = [len(arms[n][0]) for n in same]
        L.append(f'  {"identical-arm length, median":38s} {median(sl):>8.0f} nt')

    if diff and a.examples:
        L.append('')
        L.append(f'examples (first {a.examples}):')
        for n in sorted(diff)[:a.examples]:
            our_seq, off = arms[n]
            d = pub[n]
            L.append(f'  {n}')
            L.append(f'    read {len(d["seq"]):>3} nt | published guide {d["g0"]}-{d["g1"]}'
                     f' target {d["t0"]}-{d["t1"]} ({d["t1"] - d["t0"]} nt)')
            L.append(f'    ours: offset {off} length {len(our_seq)}')
            L.append(f'      theirs {d["seq"][d["t0"]:d["t1"]][:60]}')
            L.append(f'      ours   {our_seq[:60]}')

    text = '\n'.join(L)
    print(text)
    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        with open(a.out, 'w') as fh:
            fh.write(text + '\n')
        print(f'\nwrote {a.out}')


if __name__ == '__main__':
    main()
