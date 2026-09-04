#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Post-process the sno-chimeras output for the AluACA question.

The pipeline emits one `<uid>.snoRNA.<target>.chimeras.csv` per target (rRNA, snRNA,
tRNA and the genome). This script pools them and answers three things the pipeline
itself does not:

  1. Which chimeras are guided by an **AluACA** rather than a canonical snoRNA.
     In the merged source catalogue AluACA records are exactly those whose FASTA id
     ends `.id3xxx` (ids 3001-3765, numbered above snoRNA.txt.fa's highest id so the
     two sets cannot collide). A read can have several equal-scoring guides, so each
     chimera is classed AluACA / snoRNA / ambiguous rather than assigned outright.

  2. What the genomic arm actually hit -- gene, gene type and feature -- by
     intersecting the hg38 arm against a GENCODE annotation. `protein_coding` exonic
     hits are the AluACA-mRNA chimeras.

  3. Which chimeras are suspect. Two flags, kept as columns rather than applied as a
     filter, so the cost of each is visible:
       * `target_in_source_locus` - the genomic arm lands on a guide locus on the same
         strand, i.e. the "chimera" is likely one contiguous transcript.
       * `target_in_repeat` - the genomic arm lands in an annotated repeat. Alu is the
         single largest repeat family and every AluACA is Alu-derived, so an
         Alu-to-Alu pairing is the dominant false-positive mode here.
"""

import argparse
import glob
import gzip
import os
import subprocess
import sys
from collections import defaultdict

import pandas as pd


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def load_alu_names(path):
    """The exact set of AluACA record names, read from the union FASTA.

    Do not infer this from the identifier. The union records are numbered id3001-id3765,
    but 92 names in snoRNA.txt.fa (ACA64.id366, SCARNA17.id347, SNORA73B.id384, ...) also
    begin "id3", so a prefix test silently reclassifies them as AluACA. All 765 union
    names appear verbatim in the merged catalogue and none collides with a snoRNA name,
    so exact membership is both available and unambiguous.
    """
    if not path or not os.path.exists(path):
        sys.exit(f'AluACA FASTA {path!r} not found; it is needed to class guides.')
    names = {l[1:].strip() for l in open(path) if l.startswith('>')}
    if not names:
        sys.exit(f'No FASTA headers in {path}.')
    return names


def guide_class(reference, alu_names):
    """Class a '|'-joined list of guide names as AluACA / snoRNA / ambiguous."""
    # bowtie2 -a can report the same reference more than once for one read, so the
    # joined string may repeat a name; dedupe before deciding.
    names = set(str(reference).split('|'))
    alu = names & alu_names
    if len(alu) == len(names):
        return 'AluACA'
    if not alu:
        return 'snoRNA'
    return 'ambiguous'


def load_chimeras(outdir, uid, stag, tags):
    frames = []
    for tag in tags:
        path = os.path.join(outdir, f'{uid}.{stag}.{tag}.chimeras.csv')
        if not os.path.exists(path):
            print(f'  ! {os.path.basename(path)} missing, skipped', file=sys.stderr)
            continue
        df = pd.read_csv(path)
        if df.empty:
            continue
        # Per-target column names carry the tag; normalise so the pool is one schema.
        ren = {}
        for c in df.columns:
            if c.startswith(f'map_to_{tag}_'):
                ren[c] = c.replace(f'map_to_{tag}_', 'map_to_target_')
            elif c == f'reference_{tag}':
                ren[c] = 'reference_target'
        df = df.rename(columns=ren)
        df['target_class'] = tag
        frames.append(df)
        print(f'  {tag:6s} {df.shape[0]:>8,} chimeric reads')
    if not frames:
        sys.exit('No chimeras CSV found; did the pipeline finish?')
    return pd.concat(frames, ignore_index=True, sort=False)


def bed_annotate(df, gtag, gtf, source_bed, rmsk, bedtools, workdir):
    """Annotate the genomic arm of the hg38 chimeras via bedtools."""
    g = df[df.target_class == gtag].copy()
    if g.empty:
        print('  no genomic chimeras to annotate')
        for c in ['gene_name', 'gene_type', 'feature', 'target_in_source_locus', 'target_in_repeat']:
            df[c] = pd.NA
        return df

    os.makedirs(workdir, exist_ok=True)
    arm = os.path.join(workdir, 'target_arm.bed')
    g['chim_idx'] = g.index
    with open(arm, 'w') as o:
        for r in g.itertuples():
            # ref_start was made 1-based by identify_chimeric_read; back to BED.
            start = int(r.map_to_target_ref_start) - 1
            stop = int(r.map_to_target_ref_stop)
            if stop <= start:
                continue
            o.write(f'{r.reference_target}\t{start}\t{stop}\t{r.chim_idx}\t0\t{r.map_to_target_strand}\n')
    subprocess.run(f'LC_ALL=C sort -k1,1 -k2,2n {arm} -o {arm}', shell=True, check=True)

    def intersect(b_file, extra=''):
        """Return {chimera index: set(labels)} for arm x b_file."""
        cmd = f'{bedtools} intersect -a {arm} -b {b_file} -wa -wb {extra}'
        out = defaultdict(set)
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if p.returncode:
            print(f'  ! bedtools failed: {p.stderr.strip()[:200]}', file=sys.stderr)
            return out
        for line in p.stdout.splitlines():
            f = line.split('\t')
            out[int(f[3])].add(f[9] if len(f) > 9 else '1')
        return out

    # --- genes from the GENCODE GTF -----------------------------------------
    genes_bed = os.path.join(workdir, 'genes.bed')
    exons_bed = os.path.join(workdir, 'exons.bed')
    if not os.path.exists(genes_bed):
        print('  flattening GENCODE annotation ...')
        with _open(gtf) as f, open(genes_bed, 'w') as og, open(exons_bed, 'w') as oe:
            for line in f:
                if line.startswith('#'):
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 9 or c[2] not in ('gene', 'exon'):
                    continue
                attr = c[8]
                def get(k):
                    i = attr.find(k + ' "')
                    if i < 0:
                        return '.'
                    i += len(k) + 2
                    return attr[i:attr.find('"', i)]
                row = f'{c[0]}\t{int(c[3])-1}\t{c[4]}\t{get("gene_name")}\t{get("gene_type")}\t{c[6]}\n'
                (og if c[2] == 'gene' else oe).write(row)
        for b in (genes_bed, exons_bed):
            subprocess.run(f'LC_ALL=C sort -k1,1 -k2,2n {b} -o {b}', shell=True, check=True)

    print('  intersecting genomic arm with genes / exons ...')
    gene_hits = intersect(genes_bed, '-s')
    type_hits = defaultdict(set)
    p = subprocess.run(f'{bedtools} intersect -a {arm} -b {genes_bed} -wa -wb -s',
                       shell=True, capture_output=True, text=True)
    for line in p.stdout.splitlines():
        f = line.split('\t')
        type_hits[int(f[3])].add(f[10])
    exon_hits = intersect(exons_bed, '-s')

    ann = {}
    for i in g['chim_idx']:
        names = sorted(gene_hits.get(i, []))
        types = sorted(type_hits.get(i, []))
        ann[i] = (
            '|'.join(names) if names else 'intergenic',
            '|'.join(types) if types else 'intergenic',
            'exonic' if i in exon_hits else ('intronic' if names else 'intergenic'),
        )
    df.loc[g['chim_idx'], 'gene_name'] = [ann[i][0] for i in g['chim_idx']]
    df.loc[g['chim_idx'], 'gene_type'] = [ann[i][1] for i in g['chim_idx']]
    df.loc[g['chim_idx'], 'feature'] = [ann[i][2] for i in g['chim_idx']]

    # --- suspect flags -------------------------------------------------------
    if source_bed and os.path.exists(source_bed):
        print('  flagging arms that land on a guide locus ...')
        sb = os.path.join(workdir, 'source.sorted.bed')
        subprocess.run(f'LC_ALL=C sort -k1,1 -k2,2n {source_bed} -o {sb}', shell=True, check=True)
        hits = intersect(sb, '-s -f 0.25')
        df.loc[g['chim_idx'], 'target_in_source_locus'] = [i in hits for i in g['chim_idx']]
    if rmsk and os.path.exists(rmsk):
        print('  flagging arms that land in a repeat ...')
        hits = intersect(rmsk, '-f 0.5')
        df.loc[g['chim_idx'], 'target_in_repeat'] = [i in hits for i in g['chim_idx']]
    return df


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--outdir', required=True, help='Pipeline output directory.')
    p.add_argument('--uid', required=True, help='Sample identifier.')
    p.add_argument('--stag', default='snoRNA', help='Source RNA tag, default: %(default)s.')
    p.add_argument('--tags', default='rRNA,snRNA,tRNA,hg38',
                   help='Comma separated target tags, default: %(default)s.')
    p.add_argument('--gtag', default='hg38', help='Genome tag, default: %(default)s.')
    p.add_argument('--gtf', default='/home/rafail/Downloads/hg38/gencode.v47.primary_assembly.annotation.gtf.gz')
    p.add_argument('--alu-fasta', default='data/AluACA_union_nr.fasta',
                   help='FASTA whose headers name the AluACA records, default: %(default)s.')
    p.add_argument('--source-bed', default='', help='BED of guide loci, for the false-chimera flag.')
    p.add_argument('--rmsk', default='', help='BED of repeats, for the Alu-to-Alu flag.')
    p.add_argument('--bedtools', default='bedtools')
    p.add_argument('--out', required=True, help='Output TSV of annotated chimeras.')
    args = p.parse_args()

    tags = [t for t in args.tags.split(',') if t]
    print(f'Loading chimeras for {args.uid}:')
    df = load_chimeras(args.outdir, args.uid, args.stag, tags)
    alu_names = load_alu_names(args.alu_fasta)
    print(f'  {len(alu_names)} AluACA record names loaded from {args.alu_fasta}')
    df['guide_class'] = df[f'reference_{args.stag}'].map(lambda r: guide_class(r, alu_names))
    # Collapse the repeated-name artefact so guide tallies count references, not alignments.
    df['guide_names'] = df[f'reference_{args.stag}'].map(
        lambda r: '|'.join(sorted(set(str(r).split('|')))))
    print(f'  pooled {df.shape[0]:,} chimeric reads')

    df = bed_annotate(df, args.gtag, args.gtf, args.source_bed, args.rmsk,
                      args.bedtools, os.path.join(args.outdir, 'annotate_work'))
    df.to_csv(args.out, sep='\t', index=False)
    print(f'\nWrote {args.out}')

    print('\n=== chimeras by guide class and target ===')
    print(pd.crosstab(df['guide_class'], df['target_class'], margins=True).to_string())

    alu = df[df.guide_class == 'AluACA']
    if not alu.empty:
        print(f'\n=== AluACA-guided chimeras: {alu.shape[0]:,} ===')
        g = alu[alu.target_class == args.gtag]
        if not g.empty and 'gene_type' in g:
            print('\ntop genomic target biotypes:')
            print(g['gene_type'].value_counts().head(10).to_string())
            mrna = g[g.gene_type.astype(str).str.contains('protein_coding', na=False)]
            print(f'\nAluACA-mRNA chimeras (protein_coding genomic arm): {mrna.shape[0]:,}')
            if not mrna.empty:
                print('  by feature:', dict(mrna['feature'].value_counts()))
                print('\n  top mRNA targets:')
                print(mrna['gene_name'].value_counts().head(15).to_string())
        print('\ntop AluACA guides:')
        print(alu['guide_names'].value_counts().head(15).to_string())
        for flag in ('target_in_source_locus', 'target_in_repeat'):
            if flag in alu.columns:
                n = alu[flag].fillna(False).astype(bool).sum()
                print(f'\nflagged {flag}: {n:,} of {alu.shape[0]:,} AluACA chimeras')


if __name__ == '__main__':
    main()
