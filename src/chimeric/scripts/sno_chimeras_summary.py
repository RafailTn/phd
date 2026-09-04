#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Summarise a sno-chimeras run into a small self-contained HTML report.

Replacement for the `sno_chimeras_summary` console script that upstream
sno-chimeras.py invokes from the VanNostrand lab's unpublished `iToolBox`
package. Only the command-line contract is reproduced, not the original
report's layout -- the numbers come from the same two files the pipeline
already writes:

  * the ``--count`` TSV, appended to by ``identify_chimeric_read`` as
    ``<uid>.<stag>.<tag> <before_dedup> <after_dedup>`` (space separated), and
  * ``individual.chimeras.count.tsv``, written by ``parse_individual_chimeras``
    as a CSV with an ``uID,snoRNA,target,chimeras_count`` header repeated once
    per append.

Both are read defensively: the pipeline appends to them from several tasks, so
duplicate headers and partially written lines are expected rather than fatal.
"""

import argparse
import html
import os
from collections import defaultdict

from seqflow import logger


def _read_counts(path):
    """Parse the whitespace separated count file into {(uid, tag): (n1, n2)}."""
    counts = {}
    if not path or not os.path.exists(path):
        logger.warning(f'Count file {path} not found; the report will omit the chimera table.')
        return counts
    with open(path) as f:
        for line in f:
            fields = line.split()
            if len(fields) != 3:
                continue
            name, n1, n2 = fields
            try:
                n1, n2 = int(n1), int(n2)
            except ValueError:
                continue
            # name is '<uid>.<stag>.<tag>'; uid may itself contain dots.
            uid, _, tag = name.rsplit('.', 2)
            counts[(uid, tag)] = (n1, n2)
    return counts


def _read_individual(path):
    """Parse individual.chimeras.count.tsv into {(uid, tag): [(rna, count), ...]}."""
    individual = defaultdict(list)
    if not os.path.exists(path):
        return individual
    with open(path) as f:
        for line in f:
            fields = [x.strip() for x in line.split(',')]
            if len(fields) != 4 or fields[0] == 'uID':
                continue
            uid, rna, tag, count = fields
            try:
                count = int(count)
            except ValueError:
                continue
            individual[(uid, tag)].append((rna, count))
    for key in individual:
        individual[key].sort(key=lambda x: x[1], reverse=True)
    return individual


def _table(rows, headers):
    head = ''.join(f'<th>{html.escape(str(h))}</th>' for h in headers)
    body = ''.join(
        '<tr>' + ''.join(f'<td>{html.escape(str(c))}</td>' for c in row) + '</tr>'
        for row in rows
    )
    return f'<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


CSS = """
body{font-family:system-ui,-apple-system,Segoe UI,Helvetica,Arial,sans-serif;
     margin:2rem auto;max-width:60rem;line-height:1.5;color:#1a1a1a}
h1{font-size:1.4rem} h2{font-size:1.1rem;margin-top:2rem}
table{border-collapse:collapse;margin:.75rem 0;font-size:.9rem}
th,td{border:1px solid #d0d0d0;padding:.3rem .6rem;text-align:right}
th:first-child,td:first-child{text-align:left}
thead{background:#f4f4f4}
p.note{color:#555;font-size:.85rem}
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stag', required=True, help='Tag for the source RNA.')
    parser.add_argument('--tags', required=True, help='Comma separated target tags (including the genome tag).')
    parser.add_argument('--ids', required=True, help='Comma separated sample identifiers.')
    parser.add_argument('--output', required=True, help='Path to the output HTML report.')
    parser.add_argument('--count', required=True, help='Path to the chimeras count file.')
    parser.add_argument('--gtag', required=True, help='Tag for the genome.')
    args = parser.parse_args()

    tags = [t for t in args.tags.split(',') if t]
    ids = [i for i in args.ids.split(',') if i]
    counts = _read_counts(args.count)
    individual = _read_individual('individual.chimeras.count.tsv')

    sections = [f'<h1>{html.escape(args.stag)} chimeras summary</h1>']
    sections.append(
        f'<p class="note">Source RNA: <code>{html.escape(args.stag)}</code> &middot; '
        f'targets: <code>{html.escape(", ".join(tags))}</code> &middot; '
        f'genome: <code>{html.escape(args.gtag)}</code></p>'
    )

    rows = []
    for uid in ids:
        for tag in tags:
            n1, n2 = counts.get((uid, tag), (0, 0))
            rows.append([uid, tag, f'{n1:,}', f'{n2:,}'])
        total1 = sum(counts.get((uid, t), (0, 0))[0] for t in tags)
        total2 = sum(counts.get((uid, t), (0, 0))[1] for t in tags)
        rows.append([uid, 'all targets', f'{total1:,}', f'{total2:,}'])
    sections.append('<h2>Chimeric reads per target</h2>')
    sections.append(_table(rows, ['sample', 'target', 'chimeras', 'after de-duplication']))
    sections.append(
        '<p class="note">"chimeras" counts distinct read names passing the gap '
        '(|gap| &le; 4 nt) and target length (&ge; 16 nt) criteria; the second column '
        'removes UMI/position duplicates.</p>'
    )

    if individual:
        sections.append(f'<h2>Individual {html.escape(args.stag)} genes with the most chimeras</h2>')
        for uid in ids:
            for tag in tags:
                items = individual.get((uid, tag), [])
                if not items:
                    continue
                sections.append(f'<h3>{html.escape(uid)} &mdash; {html.escape(tag)}</h3>')
                sections.append(_table([[rna, f'{n:,}'] for rna, n in items[:50]],
                                       [args.stag, 'chimeras']))

    with open(args.output, 'w') as o:
        o.write(f'<!doctype html><meta charset="utf-8"><title>{html.escape(args.stag)} chimeras</title>'
                f'<style>{CSS}</style>' + ''.join(sections))
    logger.info(f'Summary report saved to {args.output}.')


if __name__ == '__main__':
    main()
