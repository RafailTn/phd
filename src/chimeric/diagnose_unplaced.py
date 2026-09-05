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

    ours_named = csv_names(os.path.join(a.outdir, f'{a.uid}.{a.stag}.{a.gtag}.chimeras.csv'))
    missed = set(pub) - ours_named
    shared = set(pub) & ours_named

    arms_path = os.path.join(a.outdir, f'{a.uid}.mask.map.to.{a.stag}.true.target.fasta')
    raw = read_fasta(arms_path)
    if not raw:
        sys.exit(f'no target arms at {arms_path} (the run needs --keep)')
    arms = {}
    for header, seq in raw.items():
        base, off = split_offset(header)
        arms[base] = (seq, off)

    def profile(names):
        """How our extracted arm relates to the region the publication aligned.

        Run over BOTH groups on purpose: a discrepancy among the unplaced reads
        only means something if the reads that placed fine do not show it too.
        """
        st = dict(n=0, exact=0, extends=0, other=0, dlen=[], dstart=[])
        for n in names:
            if n not in arms:
                continue
            our_seq, off = arms[n]
            d = pub[n]
            theirs = d['seq'][d['t0']:d['t1']]
            st['n'] += 1
            st['dlen'].append(len(our_seq) - len(theirs))
            if off is not None:
                st['dstart'].append(off - d['t0'])
            if our_seq == theirs:
                st['exact'] += 1
            elif theirs and our_seq.endswith(theirs):
                st['extends'] += 1      # same 3' end, extra bases on the 5' side
            else:
                st['other'] += 1
        return st

    m, sh = profile(missed), profile(shared)

    L = []
    L.append(f'# our target arm vs the region the publication aligned  ({a.outdir})')
    L.append('')
    L.append(f'{"published chimeras":38s} {len(pub):>8,}')
    L.append(f'{"not called by this arm (missed)":38s} {len(missed):>8,}')
    L.append(f'{"called by this arm (shared)":38s} {len(shared):>8,}')
    L.append('')
    L.append(f'{"":38s} {"missed":>10} {"shared":>10}')
    L.append(f'{"arms available to compare":38s} {m["n"]:>10,} {sh["n"]:>10,}')
    L.append(f'{"our arm == their aligned region":38s} {pct(m["exact"], m["n"]):>10} {pct(sh["exact"], sh["n"]):>10}')
    L.append(f'{"our arm extends it at the 5 prime":38s} {pct(m["extends"], m["n"]):>10} {pct(sh["extends"], sh["n"]):>10}')
    L.append(f'{"genuinely different sequence":38s} {pct(m["other"], m["n"]):>10} {pct(sh["other"], sh["n"]):>10}')
    for key, lbl in (('dlen', 'length difference, median'),
                     ('dstart', 'start offset difference, median')):
        mv = f'{median(m[key]):.0f}' if m[key] else 'n/a'
        sv = f'{median(sh[key]):.0f}' if sh[key] else 'n/a'
        L.append(f'{lbl:38s} {mv:>10} {sv:>10}')
    L.append('')
    if m['n'] and sh['n']:
        if abs(m['exact'] / m['n'] - sh['exact'] / sh['n']) < 0.05:
            L.append('The two groups look alike, so how the arm is cut is NOT what separates')
            L.append('placed from unplaced. The difference is in the alignment itself:')
            L.append('index content or STAR parameters.')
        else:
            L.append('The groups differ, so the split point is doing the separating.')
            L.append('Look upstream at the bowtie2 guide assignment, not at STAR.')

    if a.examples:
        for lbl, names in (('missed', missed), ('shared', shared)):
            L.append('')
            L.append(f'examples, {lbl} (first {a.examples}):')
            for n in sorted(x for x in names if x in arms)[:a.examples]:
                our_seq, off = arms[n]
                d = pub[n]
                theirs = d['seq'][d['t0']:d['t1']]
                L.append(f'  {n}')
                L.append(f'    read {len(d["seq"]):>3} nt | their target read {d["t0"]}-{d["t1"]}'
                         f' | our arm offset {off} len {len(our_seq)}')
                L.append(f'      theirs {theirs[:58]}')
                L.append(f'      ours   {our_seq[:58]}')

    text = '\n'.join(L)
    print(text)
    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        with open(a.out, 'w') as fh:
            fh.write(text + '\n')
        print(f'\nwrote {a.out}')


if __name__ == '__main__':
    main()
