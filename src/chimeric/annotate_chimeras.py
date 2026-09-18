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

  4. Which chimeras are pipeline artefacts. The input control is assumed not to go
     through the chimeric ligation -- the protocol does not say, see the README -- so
     every chimera called there is treated as a false positive, and these flags are
     what the input calls turned out to be:
       * `guide_low_complexity` - the guide arm's DUST score is >= 2: it is a simple
         repeat (poly(A), (GA)n, (TG)n ...) that matches too many sequences to say which
         RNA it came from. AluACA records carry the A-rich Alu tail, so an mRNA 3' end
         plus its poly(A) tail is otherwise called an AluACA chimera.
       * `genome_contiguity` - the whole read aligned against the entire genome with
         STAR: `contiguous` if one alignment anywhere covers both arms, `too many loci`
         if it maps to more places than can be listed, else `no`. Catches multi-copy
         transcripts whose "target" was placed at the wrong copy, which the local
         checks below cannot, since they only look around the reported locus.
       * `read_contiguous` - one local alignment of the whole read to the genome around
         the target arm also covers the guide arm: the read is a single transcript,
         and no ligation is needed to explain it.
       * `guide_near_target` - the guide arm occurs within 2 kb of the target arm, same
         orientation (includes `read_contiguous`); catches spliced or edited reads.
         Weak for AluACA, since an Alu lies within 2 kb of most loci by chance.
     The last two need `target_arm_aligned`, which checks the read really does align
     at the reported locus; where it does not, they are left empty.
"""

import argparse
import glob
import gzip
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

import pandas as pd

# paths.py lives one level up, shared with the analysis scripts; the repo is a
# collection of scripts rather than an installed package, so put src/ on the path.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import find_input, find_tool, out_path, proj


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
    with open(path) as fh:
        names = {l[1:].strip() for l in fh if l.startswith('>')}
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

    # Create the annotation columns up front, as object dtype holding NaN. A
    # partial .loc assignment into a column that does not exist yet makes pandas
    # build it as float64 NaN and then write strings or booleans into it, which
    # is the "Setting an item of incompatible dtype" FutureWarning -- and it
    # prints the entire value list, which for a real run is tens of thousands of
    # booleans. Rows outside the genomic arm keep NaN either way, so this only
    # changes how the column is built, not what ends up in it.
    for c in ('gene_name', 'gene_type', 'feature',
              'target_in_source_locus', 'target_in_repeat'):
        if c not in df.columns:
            df[c] = pd.Series(float('nan'), index=df.index, dtype=object)

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


# --- artefact flags -----------------------------------------------------------
# Thresholds, tuned against the input library, where every call is a false positive.
# Guide-arm coverage by the contiguity alignment is bimodal (<20% or >80%), so the
# calls do not hinge on the exact cut-offs.
DUST_MAX = 2.0      # guide arms scoring at or above this are simple repeats; costs the
                    # snoRNA control ~1.6% of its calls, removes all (GA)n / poly(A) arms
MIN_COV = 0.8       # share of an arm inside the local alignment
MIN_IDENT = 0.9     # identity over the aligned columns
FLANK = 30          # genome beyond guide-arm length either side of the target arm
NEAR = 2000         # window for guide_near_target
_RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def guide_arm(r):
    return r.sequence[int(r.map_to_snoRNA_read_start) - 1:int(r.map_to_snoRNA_read_stop)]


def dust(seq):
    """DUST score (the BLAST low-complexity measure) of a sequence.

    Counts every overlapping 3-mer and sums c*(c-1)/2 over them, normalised by the
    number of 3-mers minus one. A simple repeat reuses a few 3-mers and scores high
    (poly(A) ~11, (GA)n ~4); ordinary sequence scores well below 1."""
    n = len(seq) - 2
    if n < 2:
        return float('nan')
    c = Counter(seq[i:i + 3] for i in range(n))
    return sum(v * (v - 1) / 2 for v in c.values()) / (n - 1)


def low_complexity(seq):
    d = dust(seq)
    return d == d and d >= DUST_MAX


def _init_worker(genome_fa):
    global _FA, _AL
    import pysam
    from Bio.Align import PairwiseAligner
    _FA = pysam.FastaFile(genome_fa)
    _AL = PairwiseAligner(mode='local', match_score=2, mismatch_score=-3,
                          open_gap_score=-5, extend_gap_score=-2)


def _window(chrom, start, stop, strand, pad):
    """Genome around [start, stop) (0-based), in read orientation."""
    w = _FA.fetch(chrom, max(0, start - pad), stop + pad).upper()
    return w.translate(_RC)[::-1] if strand == '-' else w


def _align(window, query):
    """Read-side aligned intervals and identity of the best local alignment."""
    a = _AL.align(window, query)[0]
    blocks = [tuple(b) for b in a.aligned[1]]
    same = cols = 0
    for (ws, we), (qs, _) in zip(a.aligned[0], a.aligned[1]):
        cols += we - ws
        same += sum(window[ws + k] == query[qs + k] for k in range(we - ws))
    return blocks, same / max(1, cols)


def _cov(blocks, a, b):
    return sum(max(0, min(e, b) - max(s, a)) for s, e in blocks) / max(1, b - a)


def _flag_one(r):
    """(target_arm_aligned, read_contiguous, guide_near_target) for one genomic arm.

    Aligns the whole read rather than trusting the pipeline's per-arm read boundaries,
    which are not exact at the junction: position-exact comparison validated only 57%
    of target arms, this validates ~92%."""
    seq, gs, ge, ts, te, chrom, start, stop, strand = r
    try:
        w = _window(chrom, start - 1, stop, strand, ge - gs + 1 + FLANK)
    except (KeyError, ValueError):
        return None, None, None
    blocks, ident = _align(w, seq)
    if _cov(blocks, ts - 1, te) < MIN_COV:
        return False, None, None
    contiguous = _cov(blocks, gs - 1, ge) >= MIN_COV and ident >= MIN_IDENT
    if contiguous:
        return True, True, True
    guide = seq[gs - 1:ge]
    try:
        w = _window(chrom, start - 1, stop, strand, NEAR)
    except (KeyError, ValueError):
        return True, False, None
    gblocks, gident = _align(w, guide)
    near = _cov(gblocks, 0, len(guide)) >= MIN_COV and gident >= MIN_IDENT
    return True, False, near


def artefact_flags(df, gtag, genome_fa, cpus):
    """Add the pipeline-artefact columns; see point 4 of the module docstring."""
    df['guide_arm_len'] = (df.map_to_snoRNA_read_stop - df.map_to_snoRNA_read_start + 1).astype(int)
    df['guide_dust'] = [round(dust(guide_arm(r)), 3) for r in df.itertuples()]
    df['guide_low_complexity'] = df.guide_dust.ge(DUST_MAX)
    cols = ('target_arm_aligned', 'read_contiguous', 'guide_near_target')
    for c in cols:
        df[c] = pd.Series(float('nan'), index=df.index, dtype=object)
    if not genome_fa or not os.path.exists(genome_fa):
        print(f'  ! no genome FASTA at {genome_fa!r}; contiguity flags left empty', file=sys.stderr)
        return df
    if not os.path.exists(genome_fa + '.fai'):
        sys.exit(f'{genome_fa} has no .fai index; run `samtools faidx {genome_fa}` first.')
    g = df[df.target_class == gtag]
    if g.empty:
        return df
    print(f'  aligning {len(g):,} genomic chimeras back to the genome ({cpus} processes) ...')
    rows = [(r.sequence, int(r.map_to_snoRNA_read_start), int(r.map_to_snoRNA_read_stop),
             int(r.map_to_target_read_start), int(r.map_to_target_read_stop),
             r.reference_target, int(r.map_to_target_ref_start), int(r.map_to_target_ref_stop),
             r.map_to_target_strand) for r in g.itertuples()]
    with Pool(cpus, initializer=_init_worker, initargs=(genome_fa,)) as pool:
        res = pool.map(_flag_one, rows, chunksize=200)
    for c, vals in zip(cols, zip(*res)):
        df.loc[g.index, c] = list(vals)
    return df


# --- genome-wide contiguity ----------------------------------------------------
# The pipeline's genome mask aligns reads end-to-end, so a transcript with a few
# non-genomic bases at its ends, or ~95% identical to several copies, escapes it and is
# split into a "chimera". The local checks above only look around the reported target
# locus, so a multi-copy transcript whose "target" was placed at the wrong copy passes.
#
# Two stages. STAR proposes every locus the whole read aligns to (local mode, so ragged
# ends clip). STAR's own alignment cannot make the call: it soft-clips a short guide arm
# carrying a couple of mismatches, even where the read is contiguous. So each proposed
# locus, plus the pipeline's reported one, is re-aligned with the same aligner and
# thresholds as read_contiguous -- which makes this a strict superset of that check.
MULTIMAP_MAX = 500
MAX_CANDIDATES = 20
_CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')


def _ref_span(cigar):
    return sum(int(n) for n, op in _CIGAR.findall(cigar) if op in 'MDN=X')


def _contig_anywhere(r):
    """(contiguous?, loci STAR proposed) for one call, over all candidate loci."""
    seq, gs, ge, ts, te, cands, n_loci = r
    pad = len(seq) + FLANK
    for chrom, start, stop, strand in cands:
        try:
            w = _window(chrom, start, stop, strand, pad)
        except (KeyError, ValueError):
            continue
        blocks, ident = _align(w, seq)
        if _cov(blocks, gs - 1, ge) >= MIN_COV and _cov(blocks, ts - 1, te) >= MIN_COV \
                and ident >= MIN_IDENT:
            return 'contiguous', n_loci
    return 'no', n_loci


def genome_contiguity(df, star, index, genome_fa, cpus, workdir):
    """Add `genome_contiguity` and `genome_loci`; see the comment block above."""
    df['genome_contiguity'] = pd.Series(float('nan'), index=df.index, dtype=object)
    df['genome_loci'] = pd.Series(float('nan'), index=df.index, dtype=object)
    if not index or not os.path.exists(os.path.join(index, 'SA')):
        print(f'  ! no STAR index at {index!r}; genome-wide contiguity left empty', file=sys.stderr)
        return df
    if not genome_fa or not os.path.exists(genome_fa + '.fai'):
        print(f'  ! genome-wide contiguity needs the genome FASTA with .fai too; left empty',
              file=sys.stderr)
        return df
    os.makedirs(workdir, exist_ok=True)
    uniq = {s: i for i, s in enumerate(pd.unique(df.sequence.astype(str)))}
    fa = os.path.join(workdir, 'reads.fa')
    with open(fa, 'w') as o:
        o.writelines(f'>{i}\n{s}\n' for s, i in uniq.items())
    print(f'  proposing loci for {len(uniq):,} distinct reads with STAR ...')
    cmd = [star, '--runMode', 'alignReads', '--genomeDir', index, '--readFilesIn', fa,
           '--runThreadN', str(cpus), '--outFileNamePrefix', os.path.join(workdir, ''),
           '--outSAMtype', 'SAM', '--outSAMunmapped', 'Within', '--outSAMattributes', 'NH',
           '--outSAMmode', 'NoQS',
           '--alignEndsType', 'Local', '--outFilterScoreMinOverLread', '0',
           '--outFilterMatchNminOverLread', '0', '--outFilterMatchNmin', '16',
           '--outFilterMismatchNoverLmax', '0.1',
           '--outFilterMultimapNmax', str(MULTIMAP_MAX),
           '--winAnchorMultimapNmax', str(2 * MULTIMAP_MAX),
           '--outSAMmultNmax', str(MAX_CANDIDATES),
           # a spliced read is still one transcript, but only across canonical introns of
           # plausible size -- otherwise two loci on one chromosome would pass as "spliced"
           '--alignIntronMax', '100000', '--outFilterIntronMotifs', 'RemoveNoncanonical']
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode:
        sys.exit(f'STAR failed ({p.returncode}): {p.stderr.strip()[-500:]}')

    cands = defaultdict(list)   # read id -> [(chrom, start0, stop, strand)]
    loci, too_many = {}, set()
    with open(os.path.join(workdir, 'Aligned.out.sam')) as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            f = line.rstrip('\n').split('\t')
            rid, flag = int(f[0]), int(f[1])
            tags = {t.split(':', 2)[0]: t.split(':', 2)[2] for t in f[11:]}
            if flag & 4:
                if tags.get('uT') == '3':
                    too_many.add(rid)
                continue
            start = int(f[3]) - 1
            cands[rid].append((f[2], start, start + _ref_span(f[5]), '-' if flag & 16 else '+'))
            loci[rid] = int(tags.get('NH', 1))

    rows, idx = [], []
    for r in df.itertuples():
        rid = uniq[str(r.sequence)]
        c = list(cands.get(rid, []))
        if r.target_class == df.attrs.get('gtag'):
            # the pipeline's own placement, so this can only add to read_contiguous
            c.append((r.reference_target, int(r.map_to_target_ref_start) - 1,
                      int(r.map_to_target_ref_stop), r.map_to_target_strand))
        if rid in too_many and not c:
            df.at[r.Index, 'genome_contiguity'] = 'too many loci'
            continue
        rows.append((str(r.sequence), int(r.map_to_snoRNA_read_start), int(r.map_to_snoRNA_read_stop),
                     int(r.map_to_target_read_start), int(r.map_to_target_read_stop),
                     c, loci.get(rid, float('nan'))))
        idx.append(r.Index)
    print(f'  re-aligning at {sum(len(x[5]) for x in rows):,} candidate loci ({cpus} processes) ...')
    with Pool(cpus, initializer=_init_worker, initargs=(genome_fa,)) as pool:
        res = pool.map(_contig_anywhere, rows, chunksize=100)
    df.loc[idx, 'genome_contiguity'] = [x[0] for x in res]
    df.loc[idx, 'genome_loci'] = [x[1] for x in res]
    # Never weaker than the local check: the wider window here can let a different local
    # alignment win (~0.3% of locally contiguous calls), so fold those in explicitly.
    if 'read_contiguous' in df:
        df.loc[df.read_contiguous.eq(True), 'genome_contiguity'] = 'contiguous'
    return df


def default_gtf(gtag):
    """The annotation fetch_refs.sh downloads, under ref/chimeric/<build>/.

    A bare relative default only works when the cwd happens to be the project
    root, and an absolute one only works on the machine it was written on; this
    is derived from the repo instead.
    """
    ref = os.environ.get('REF') or os.path.join(proj(), 'ref', 'chimeric')
    names = {'hg38': 'gencode.v47.primary_assembly.annotation.gtf.gz',
             'hg19': 'gencode.v47lift37.annotation.gtf.gz'}
    return os.path.join(ref, gtag, names.get(gtag, names['hg38']))


def default_genome_fa(gtag):
    """The primary assembly fetch_refs.sh downloads, next to the GTF."""
    names = {'hg38': 'GRCh38.primary_assembly.genome.fa',
             'hg19': 'GRCh37.primary_assembly.genome.fa'}
    return os.path.join(os.path.dirname(default_gtf(gtag)), names.get(gtag, names['hg38']))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--outdir', help='Pipeline output directory; required unless --annotated.')
    p.add_argument('--uid', help='Sample identifier; required unless --annotated.')
    p.add_argument('--stag', default='snoRNA', help='Source RNA tag, default: %(default)s.')
    p.add_argument('--tags', default='rRNA,snRNA,tRNA,hg38',
                   help='Comma separated target tags, default: %(default)s.')
    p.add_argument('--gtag', default='hg38', help='Genome tag, default: %(default)s.')
    p.add_argument('--gtf', default=None,
                   help='GENCODE annotation GTF, for gene and exon context. '
                        'Default: $GENCODE, else the one fetch_refs.sh puts in '
                        'ref/chimeric/<gtag>/.')
    p.add_argument('--alu-fasta', default=find_input('AluACA_union_nr.fasta'),
                   help='FASTA whose headers name the AluACA records, default: %(default)s.')
    p.add_argument('--source-bed', default=os.environ.get('GUIDE_BED', ''),
                   help='BED of guide loci, for the false-chimera flag.')
    p.add_argument('--rmsk', default=os.environ.get('RMSK_BED', ''),
                   help='BED of repeats, for the Alu-to-Alu flag.')
    p.add_argument('--bedtools', default=None,
                   help='bedtools executable. Default: $BEDTOOLS, the project pixi env, then PATH.')
    p.add_argument('--genome-fa', default=None,
                   help='Genome FASTA (with .fai) for the contiguity flags. Default: $GENOME_FA, '
                        'else the one fetch_refs.sh puts in ref/chimeric/<gtag>/.')
    p.add_argument('--genome-index', default=os.environ.get('GENOME_INDEX', ''),
                   help='STAR index of the genome, for the genome-wide contiguity check. '
                        'Default: $GENOME_INDEX. Needs the RAM to load it.')
    p.add_argument('--star', default=None,
                   help='STAR executable. Default: $STAR, the project pixi env, then PATH.')
    p.add_argument('--cpus', type=int, default=int(os.environ.get('CPUS') or os.cpu_count() or 1),
                   help='Processes for the contiguity alignments, default: %(default)s.')
    p.add_argument('--annotated', default='',
                   help='Existing annotated TSV: only (re)compute the artefact flags on it, '
                        'skipping pooling and bedtools. For results whose per-target CSVs are '
                        'elsewhere.')
    p.add_argument('--out', required=True, help='Output TSV of annotated chimeras.')
    args = p.parse_args()
    if not args.annotated and not (args.outdir and args.uid):
        p.error('--outdir and --uid are required unless --annotated is given')
    # Resolved here rather than as argparse defaults: the GTF depends on --gtag,
    # and the tool lookup should not run when an explicit path was given.
    args.gtf = args.gtf or os.environ.get('GENCODE') or default_gtf(args.gtag)
    args.genome_fa = args.genome_fa or os.environ.get('GENOME_FA') or default_genome_fa(args.gtag)
    args.star = find_tool('STAR', args.star)

    if args.annotated:
        df = pd.read_csv(args.annotated, sep='\t', low_memory=False)
        print(f'Re-flagging {df.shape[0]:,} chimeras from {args.annotated}')
        df = artefact_flags(df, args.gtag, args.genome_fa, args.cpus)
        df.attrs['gtag'] = args.gtag
        df = genome_contiguity(df, args.star, args.genome_index, args.genome_fa, args.cpus,
                               os.path.abspath(args.out) + '.genome_work')
        df.to_csv(args.out, sep='\t', index=False)
        print(f'\nWrote {args.out}')
        summarise_flags(df)
        return
    args.bedtools = find_tool('bedtools', args.bedtools)

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
    df = artefact_flags(df, args.gtag, args.genome_fa, args.cpus)
    df.attrs['gtag'] = args.gtag
    df = genome_contiguity(df, args.star, args.genome_index, args.genome_fa, args.cpus,
                           os.path.join(args.outdir, 'genome_contiguity_work'))
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
                # .eq rather than .fillna(False).astype(bool): NaN.eq(True) is
                # already False, and filling an object column then letting astype
                # downcast it is deprecated in pandas 2.x.
                n = alu[flag].eq(True).sum()
                print(f'\nflagged {flag}: {n:,} of {alu.shape[0]:,} AluACA chimeras')
    summarise_flags(df)


def summarise_flags(df):
    """Share of each guide class's genomic chimeras carrying each artefact flag."""
    g = df[df.target_arm_aligned.eq(True)] if 'target_arm_aligned' in df else df.iloc[:0]
    if g.empty:
        return
    gw = g.genome_contiguity.isin(['contiguous', 'too many loci'])
    any_flag = g.guide_low_complexity.eq(True) | g.guide_near_target.eq(True) | gw
    t = pd.DataFrame({
        'aligned': g.groupby('guide_class').size(),
        'low_complexity': g.guide_low_complexity.eq(True).groupby(g.guide_class).mean(),
        'contiguous': g.read_contiguous.eq(True).groupby(g.guide_class).mean(),
        'near_target': g.guide_near_target.eq(True).groupby(g.guide_class).mean(),
        'genome-wide': gw.groupby(g.guide_class).mean(),
        'unflagged': (~any_flag).groupby(g.guide_class).mean(),
    })
    print('\n=== artefact flags, genomic chimeras whose target arm re-aligned ===')
    print(t.to_string(float_format=lambda v: f'{v:.1%}'))


if __name__ == '__main__':
    main()
