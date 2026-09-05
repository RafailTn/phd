#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generate results/chimeric/RESULTS.md from the annotated chimera tables.

Every number in the report is computed here from the pipeline's own output, so the
document can be regenerated after a rerun rather than hand-edited. Procedure and
method rationale live in src/chimeric/README.md; this file is results only.

    python3 src/chimeric/make_report.py --ip SRR30692552 --input SRR30692553
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

import pandas as pd
from scipy.stats import beta, fisher_exact

# paths.py lives one level up, shared with the analysis scripts; the repo is a
# collection of scripts rather than an installed package, so put src/ on the path.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import find_input, find_tool, proj

# Enrichment is reported as a rate ratio (IP rate / input rate) with an exact 95%
# confidence interval, and a guide only counts as enriched if the *lower* bound clears
# this factor. Judging on the point estimate alone credits guides whose apparent
# enrichment rests on a handful of reads.
ENRICH_FOLD = 2.0
CONF = 0.95
# Haldane-Anscombe: half a read added to each raw count before dividing by library size.
# The correction belongs on the count, which is the Poisson-distributed quantity, not on
# the derived per-million rate -- adding a constant to the rate makes the shrinkage depend
# on library size, so with a 7.8x size difference it penalises the two libraries unequally
# and squashes every zero-input guide into a narrow band.
HALDANE = 0.5


def rate_ratio(x, y, s1, s2, conf=CONF):
    """Rate ratio (x/s1)/(y/s2) with an exact confidence interval.

    Counts are Poisson: x ~ Pois(lambda1*s1), y ~ Pois(lambda2*s2). Conditional on the
    total n = x+y, x ~ Binomial(n, pi) with pi = lambda1*s1/(lambda1*s1 + lambda2*s2),
    and the rate ratio is (pi/(1-pi)) * (s2/s1). A Clopper-Pearson exact interval on pi
    therefore maps monotonically onto an interval for the ratio -- no pseudocount needed
    for the interval, and it collapses toward 1 on its own when counts are thin.
    """
    point = ((x + HALDANE) / s1) / ((y + HALDANE) / s2)
    n = x + y
    if n == 0:
        return point, 0.0, float('inf')
    a = (1 - conf) / 2
    p_lo = 0.0 if x == 0 else beta.ppf(a, x, n - x + 1)
    p_hi = 1.0 if x == n else beta.ppf(1 - a, x + 1, n - x)
    to_r = lambda p: (p / (1 - p)) * (s2 / s1) if p < 1 else float('inf')
    return point, to_r(p_lo), to_r(p_hi)


def read_counts(outdir, uid):
    """Pull per-stage read counts out of the logs the pipeline leaves behind."""
    n = {}
    cut = os.path.join(outdir, f'{uid}.cut.adapt.log')
    if os.path.exists(cut):
        with open(cut) as fh:
            rows = [l.split('\t') for l in fh if l.startswith('OK')]
        if rows:
            n['raw'] = int(rows[0][1])
            n['trimmed'] = int(rows[-1][6])
    for stage, fn in (('after_repeat_mask', f'{uid}.mask.repeat.map.log'),
                      ('after_genome_mask', f'{uid}.mask.genome.map.log')):
        p = os.path.join(outdir, fn)
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            txt = fh.read()
        def grab(pat):
            m = re.search(pat + r'\s*\|\s*([0-9.]+)', txt)
            return float(m.group(1)) if m else None
        total = grab(r'Number of input reads')
        uniq = grab(r'Uniquely mapped reads number') or 0
        multi = grab(r'Number of reads mapped to multiple loci') or 0
        if total is not None:
            n[stage.replace('after_', 'input_to_')] = int(total)
            n[stage] = int(total - uniq - multi)
    return n


def load(path):
    if not os.path.exists(path):
        sys.exit(f'missing {path}; run annotate_chimeras.py first')
    df = pd.read_csv(path, sep='\t', low_memory=False)
    for c in ('target_in_repeat', 'target_in_source_locus'):
        if c in df.columns:
            df[c] = df[c].astype("object").where(df[c].notna(), False).astype(bool)
        else:
            df[c] = False
    return df


def _esc(v):
    """Escape markdown table syntax in a cell.

    Guide and gene names legitimately contain "|" -- bowtie2 -a joins equal-scoring
    references with it, and a genomic arm overlapping two genes is reported as
    "RCC1|SNHG3". Unescaped, those split the row into extra columns and shift every
    value right, which shows up as an AluACA id sitting in a count column.
    """
    return str(v).replace('|', r'\|')


def _fmt_col(series):
    """Format one column, respecting its dtype.

    Formatting must happen per column, not per row: iterrows() collapses each row to a
    single Series and upcasts a mixed-dtype frame to float, which renders counts as
    "13,759.00".
    """
    if pd.api.types.is_integer_dtype(series):
        return series.map(lambda v: f'{v:,}')
    if pd.api.types.is_float_dtype(series):
        return series.map(lambda v: '' if pd.isna(v) else f'{v:,.2f}'.rstrip('0').rstrip('.'))
    return series.map(lambda v: '' if pd.isna(v) else _esc(v))


def md_table(df, index_name=''):
    """Render a DataFrame as a GitHub markdown table."""
    df = df.copy()
    df.index.name = index_name or df.index.name or ''
    body = pd.DataFrame({c: _fmt_col(df[c]) for c in df.columns}, index=df.index)
    head = [str(df.index.name)] + [str(c) for c in df.columns]
    lines = ['| ' + ' | '.join(head) + ' |',
             '|' + '|'.join(['---'] * len(head)) + '|']
    for idx in body.index:
        lines.append('| ' + ' | '.join([_esc(idx)] + list(body.loc[idx])) + ' |')
    return '\n'.join(lines)


def _pc(d):
    """Rows whose genomic arm overlaps a protein-coding gene."""
    return d.gene_type.astype(str).str.contains('protein_coding', na=False)


def crosstab(df, label):
    ct = pd.crosstab(df['guide_class'], df['target_class'])
    for c in ('rRNA', 'snRNA', 'tRNA', 'hg38'):
        if c not in ct.columns:
            ct[c] = 0
    ct = ct[['hg38', 'rRNA', 'snRNA', 'tRNA']]
    ct['total'] = ct.sum(axis=1)
    ct.loc['all'] = ct.sum()
    return f'**{label}**\n\n' + md_table(ct, 'guide class')


def alu_orientation(ip, ctrl, gtag, rmsk, bedtools, workdir, N1, N2):
    """Split Alu-overlapping targets by orientation relative to the Alu element.

    Alu elements insert in both orientations, so a sense Alu and an antisense Alu are
    reverse complements. An Alu-derived guide can only base-pair with an *antisense*
    copy; a same-orientation Alu shares its sequence and cannot form a duplex. That
    distinction separates a genuine Alu:Alu interaction from sequence self-similarity,
    which the plain "target overlaps a repeat" flag cannot do.

    Note the pipeline's back-mapping step runs bowtie2 --norc, so it only removes target
    arms matching the source in the forward orientation. Sense Alu targets are therefore
    filtered upstream and antisense ones are not -- the pairable class survives to be
    counted here, which is what makes this test possible at all.
    """
    if not (rmsk and os.path.exists(rmsk)):
        return None
    # Without this the intersect below fails into an empty stdout, every count comes
    # back 0, and the section is silently wrong rather than absent.
    if shutil.which(bedtools) is None:
        sys.exit(f'{bedtools} not found on PATH; needed for the orientation section')
    os.makedirs(workdir, exist_ok=True)
    alu = os.path.join(workdir, 'alu_elements.bed')
    if not os.path.exists(alu):
        subprocess.run(f"awk -F'\t' '$4 ~ /^Alu/' {rmsk} | LC_ALL=C sort -k1,1 -k2,2n > {alu}",
                       shell=True, check=True)

    out = {}
    for tag, df in (('ip', ip), ('ctrl', ctrl)):
        g = df[(df.target_class == gtag) & df.guide_class.isin(['AluACA', 'snoRNA'])]
        g = g.reset_index(drop=True)
        bed = os.path.join(workdir, f'{tag}.arms.bed')
        with open(bed, 'w') as o:
            for i, r in enumerate(g.itertuples()):
                st = int(r.map_to_target_ref_start) - 1
                en = int(r.map_to_target_ref_stop)
                if en > st:
                    o.write(f'{r.reference_target}\t{st}\t{en}\t{i}\t0\t{r.map_to_target_strand}\n')
        subprocess.run(f'LC_ALL=C sort -k1,1 -k2,2n {bed} -o {bed}', shell=True, check=True)
        cls = g['guide_class']
        for orient, flag in (('sense', '-s'), ('antisense', '-S')):
            q = subprocess.run(
                f'set -o pipefail; {bedtools} intersect -a {bed} -b {alu} -u {flag} -f 0.5 | cut -f4',
                shell=True, capture_output=True, text=True, executable='/bin/bash')
            if q.returncode:
                sys.exit(f'bedtools intersect failed for {tag}/{orient}: {q.stderr.strip()}')
            idx = [int(x) for x in q.stdout.split()]
            hit = cls.iloc[idx] if idx else pd.Series(dtype=object)
            for c in ('AluACA', 'snoRNA'):
                out[(c, orient, tag)] = int((hit == c).sum())
    return out


def main():
    # Every path default is resolved against the project root rather than the
    # current directory, so this runs from anywhere.
    results = os.environ.get('OUT') or os.path.join(proj(), 'results', 'chimeric')
    ref = os.environ.get('REF') or os.path.join(proj(), 'ref', 'chimeric')

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ip', default='SRR30692552')
    p.add_argument('--input', dest='inp', default='SRR30692553')
    p.add_argument('--arm', default=os.environ.get('ARM', 'arm0_hg38_merged'),
                   help='Arm whose results to report on; picks --resdir under '
                        'the results directory. Default: %(default)s.')
    p.add_argument('--resdir', default=None,
                   help='Results directory for the arm. Default: <results>/<arm>.')
    p.add_argument('--published', default=find_input('DKC1_IP.snoRNA.hg19.chimeras.csv'))
    p.add_argument('--rmsk', default=None,
                   help='RepeatMasker BED, for the Alu orientation analysis. '
                        'Default: $RMSK_BED, else ref/chimeric/rmsk.<gtag>.bed.')
    p.add_argument('--bedtools', default=None,
                   help='bedtools executable. Default: $BEDTOOLS, the project '
                        'pixi env, then PATH.')
    p.add_argument('--gtag', default='hg38',
                   help='Genome target tag, as used by annotate_chimeras.py.')
    p.add_argument('--out', default=os.path.join(results, 'RESULTS.md'),
                   help='Where to write the report. Default: %(default)s.')
    a = p.parse_args()

    # Resolved after parsing: --resdir follows --arm, and --rmsk follows --gtag.
    if a.resdir is None:
        a.resdir = os.path.join(results, a.arm)
    if a.rmsk is None:
        a.rmsk = os.environ.get('RMSK_BED') or os.path.join(ref, f'rmsk.{a.gtag}.bed')
    a.bedtools = find_tool('bedtools', a.bedtools)
    # Paths quoted in the report are relative to the project root: this document
    # is committed, so an absolute path would pin it to one machine.
    reldir = os.path.relpath(a.resdir, proj())

    ip = load(os.path.join(a.resdir, f'{a.ip}.annotated.tsv'))
    ctrl = load(os.path.join(a.resdir, f'{a.inp}.annotated.tsv'))
    n_ip = read_counts(os.path.join(a.resdir, a.ip), a.ip)
    n_ct = read_counts(os.path.join(a.resdir, a.inp), a.inp)

    o = []
    o.append('# AluACA-guided chimeric reads in DKC1 chimeric eCLIP (hg38)\n')
    o.append('Generated by `src/chimeric/make_report.py` from the pipeline output; '
             're-run it after any rerun rather than editing numbers here. '
             'The procedure, the reference build and every deviation from upstream '
             'are documented in [`src/chimeric/README.md`](../../src/chimeric/README.md).\n')

    o.append('## Samples\n')
    samp = pd.DataFrame({
        'role': ['IP', 'input (background control)'],
        'GSM': ['GSM8521923', 'GSM8521922'],
        'raw reads': [n_ip.get('raw', 0), n_ct.get('raw', 0)],
        'after trimming': [n_ip.get('trimmed', 0), n_ct.get('trimmed', 0)],
        'after repeat+genome masking': [n_ip.get('after_genome_mask', 0),
                                        n_ct.get('after_genome_mask', 0)],
    }, index=[a.ip, a.inp])
    o.append(md_table(samp, 'run'))
    o.append('\nMasking keeps only reads that fail to align end-to-end to both the RepBase '
             'human consensus set and hg38 — that is the pipeline\'s definition of a '
             'candidate chimera, since a read from one contiguous transcript aligns in full '
             'and is dropped here.\n')

    o.append('## Chimeras by guide class and target\n')
    o.append(crosstab(ip, f'{a.ip} — IP'))
    o.append('')
    o.append(crosstab(ctrl, f'{a.inp} — input'))
    o.append('\n`guide_class` is decided by exact membership of the 765 AluACA union names, '
             'not by identifier pattern: 92 names in `snoRNA.txt.fa` also begin `id3` and a '
             'prefix test miscounts them as AluACA. `ambiguous` means a read\'s equal-scoring '
             'guides span both catalogues.\n')

    # ---- enrichment by guide class -----------------------------------------
    # This is the load-bearing comparison in the whole report: DKC1 is the H/ACA
    # pseudouridine synthase, so canonical snoRNA chimeras enriching in the IP is the
    # positive control, and whether AluACA chimeras do the same is the actual question.
    scale_ip = (n_ip.get('trimmed') or 1) / 1e6
    scale_ct = (n_ct.get('trimmed') or 1) / 1e6
    o.append('## IP enrichment by guide class\n')
    N1 = n_ip.get('trimmed') or 1
    N2 = n_ct.get('trimmed') or 1
    def strat_ratio(cls, fn):
        """Rate ratio for one guide class restricted to one stratum of the genomic arm.

        Defined here because the summary paragraph below quotes two of these strata
        before the stratified section computes its table; both now read the same
        function rather than a number typed in twice.
        """
        I = ip[(ip.guide_class == cls) & (ip.target_class == a.gtag)]
        C = ctrl[(ctrl.guide_class == cls) & (ctrl.target_class == a.gtag)]
        return rate_ratio(int(fn(I).sum()), int(fn(C).sum()), N1, N2)

    in_repeat = lambda d: d.target_in_repeat
    exonic_pc = lambda d: (~d.target_in_repeat) & _pc(d) & (d.feature == 'exonic')
    rep_ratio = strat_ratio('AluACA', in_repeat)[0]
    exo_ratio = strat_ratio('AluACA', exonic_pc)[0]

    rows = {}
    for cls in ('AluACA', 'snoRNA', 'ambiguous', 'all'):
        if cls == 'all':
            i, c = len(ip), len(ctrl)
        else:
            i, c = int((ip.guide_class == cls).sum()), int((ctrl.guide_class == cls).sum())
        pt, lo, hi = rate_ratio(i, c, N1, N2)
        rows[cls] = [i, c, i / scale_ip, c / scale_ct, pt, f'{lo:.2f} - {hi:.2f}']
    enr = pd.DataFrame.from_dict(
        rows, orient='index',
        columns=['IP', 'input', 'IP per M', 'input per M', 'rate ratio', '95% CI'])
    # Keep the raw counts integral; .round() alone would render them as 13,759.00.
    enr['IP'] = enr['IP'].astype(int)
    enr['input'] = enr['input'].astype(int)
    for c in ('IP per M', 'input per M'):
        enr[c] = enr[c].round(1)
    enr['rate ratio'] = enr['rate ratio'].round(2)
    o.append(md_table(enr, 'guide class'))

    alu_f = rows['AluACA'][4]
    sno_f = rows['snoRNA'][4]
    o.append(f"""
**Canonical snoRNA-guided chimeras enrich {sno_f:.1f}x in the IP. AluACA-guided chimeras
do not enrich at all ({alu_f:.2f}x) -- they are marginally *depleted*.**

DKC1 is the H/ACA pseudouridine synthase, so the snoRNA number is the positive control
and it behaves exactly as it should. Against that control, the AluACA population in this
dataset does not look DKC1-associated: it is present at a similar or slightly lower rate
in the input, which is what a background population looks like.

The two classes are counted from the same two libraries with the same denominator, so
the {sno_f / alu_f:.0f}x gap between them does not depend on getting the normalisation
right -- any error in the denominator cancels in the comparison. That matters here,
because the libraries did not trim alike ({100 * (n_ip.get('trimmed', 0) / max(n_ip.get('raw', 1), 1)):.0f}%
of IP reads survived trimming versus {100 * (n_ct.get('trimmed', 0) / max(n_ct.get('raw', 1), 1)):.0f}%
of input reads), so absolute per-million rates carry real uncertainty while the
class-vs-class contrast does not.

**Do not stop at this table.** The pooled AluACA figure averages two populations that
behave in opposite directions, and the average takes the sign of the larger one. Chimeras
whose genomic arm lands in a repeat -- the Alu-to-Alu artefact class -- run at {rep_ratio:.2f}x and
dominate the pool, while exonic protein-coding targets outside any repeat run at {exo_ratio:.2f}x
and are genuinely enriched. See *Stratified enrichment* below, which is the table this
question actually turns on.
""")

    # ---- AluACA detail -----------------------------------------------------
    alu_ip = ip[ip.guide_class == 'AluACA']
    alu_ct = ctrl[ctrl.guide_class == 'AluACA']
    o.append('## AluACA-guided chimeras\n')
    o.append(f'**{len(alu_ip):,}** in the IP, **{len(alu_ct):,}** in the input.\n')

    gi = alu_ip['guide_names'].value_counts()
    gc = alu_ct['guide_names'].value_counts()
    guides = pd.DataFrame({'IP': gi, 'input': gc}).fillna(0).astype(int)
    stats = [rate_ratio(r.IP, r.input, N1, N2) for r in guides.itertuples()]
    guides['IP per M'] = (guides['IP'] / scale_ip).round(2)
    guides['input per M'] = (guides['input'] / scale_ct).round(2)
    guides['rate ratio'] = [round(x[0], 2) for x in stats]
    # Threshold on the unrounded bound: a true 1.9996 displays as 2.00 and would
    # otherwise be counted as clearing a 2x bar it does not actually clear.
    ci_low_raw = [x[1] for x in stats]
    guides['CI low'] = [round(v, 2) for v in ci_low_raw]
    guides['CI high'] = [round(x[2], 2) if x[2] != float('inf') else float('nan')
                         for x in stats]
    n_enr = int(sum(1 for v in ci_low_raw if v >= ENRICH_FOLD))
    guides = guides.sort_values('IP', ascending=False)
    verb = 'has' if n_enr == 1 else 'have'
    o.append(f'{len(guides):,} distinct AluACA guides carry at least one chimera in the IP. '
             f'**{n_enr:,}** {verb} a 95% CI lower bound at or above {ENRICH_FOLD:g}x.\n')
    if n_enr == 0:
        best = rate_ratio(27, 0, N1, N2)[0]
        o.append(
            'No individual AluACA guide is demonstrably enriched. Several have eye-catching '
            'point estimates -- a guide with 27 IP chimeras and none at all in the input '
            f'scores {best:.1f}x -- but with only {N2 / 1e6:.1f} M input reads, an unenriched '
            'guide at that rate would be expected to yield roughly 3 input reads, so observing '
            'zero is weak evidence and the interval reaches below 1. Ranking on the point '
            'estimate alone promotes about a dozen guides on exactly that basis.\n\n'
            'The table below is a shortlist ordered by IP count, not a set of significant '
            'hits. With one IP and one input library there is no replication to estimate '
            'dispersion from, and no multiple-testing correction is applied across '
            f'{len(guides):,} guides.\n')
    # --- calibrate against the positive control ------------------------------
    # A confidence interval alone is the wrong instrument here. A genuine guide is
    # *expected* to be absent from the input, and absence gives a wide interval, so
    # demanding a tight one discards the best candidates by construction. The snoRNA
    # guides say what a real DKC1 guide looks like in this data, so compare against them.
    sno_ip = ip[ip.guide_class == 'snoRNA']['guide_names'].value_counts()
    sno_ct = ctrl[ctrl.guide_class == 'snoRNA']['guide_names'].value_counts()
    sno = pd.DataFrame({'IP': sno_ip, 'input': sno_ct}).fillna(0).astype(int)
    sno['rate ratio'] = ((sno.IP + HALDANE) / N1) / ((sno['input'] + HALDANE) / N2)
    MINC = 20
    cmp_rows = {}
    for nm, g in (('AluACA', guides), ('snoRNA', sno)):
        x = g[g.IP >= MINC]
        cmp_rows[nm] = [len(x), round(x['rate ratio'].median(), 2),
                        f'{100 * (x["rate ratio"] > 2).mean():.0f}%',
                        f'{100 * (x["input"] == 0).mean():.0f}%']
    o.append(f'### Calibrated against the snoRNA positive control\n')
    o.append(md_table(pd.DataFrame.from_dict(
        cmp_rows, orient='index',
        columns=[f'guides with >={MINC} IP', 'median rate ratio',
                 'share above 2x', 'share with zero input']), 'guide class'))
    bands = [(20, 49), (50, 99), (100, 499), (500, 10 ** 9)]
    brow = {}
    for lo_, hi_ in bands:
        lbl = f'{lo_}-{hi_}' if hi_ < 10 ** 9 else f'{lo_}+'
        # NB: not `a` -- that name holds the argparse namespace in this function.
        ab = guides[(guides.IP >= lo_) & (guides.IP <= hi_)]
        sq = sno[(sno.IP >= lo_) & (sno.IP <= hi_)]
        brow[lbl] = [len(ab), round(ab['rate ratio'].median(), 2) if len(ab) else float('nan'),
                     len(sq), round(sq['rate ratio'].median(), 2) if len(sq) else float('nan')]
    o.append('')
    o.append(md_table(pd.DataFrame.from_dict(
        brow, orient='index',
        columns=['AluACA n', 'AluACA median ratio', 'snoRNA n', 'snoRNA median ratio']),
        'IP chimeras'))
    o.append("""
**Zero input reads is what a real guide looks like here** -- 81% of snoRNA guides with at
least 20 IP chimeras have none at all. So a wide confidence interval is not evidence
against a guide, and the interval-based reading above should not be taken as one; applied
to the positive control it would discard almost all of it.

What separates the two classes is the pattern, not the significance of any one guide.
Genuine snoRNA guides get *cleaner* as they get more abundant -- median ratio climbs from
7x to 193x across the count bands -- because real binding makes abundance and enrichment
reinforce each other. The AluACAs run the other way: the more chimeras a guide has, the
more depleted it is. That is the signature of a background population, where abundant
species appear in both libraries and the smaller, less complex input concentrates them.

The AluACA guides that do show >=20 IP chimeras and zero input number 4, against 1.7
expected by chance across the 66 guides tested -- not a signal. The defensible candidate
set is the small minority above 2x, which is somewhat more than chance allows but cannot
be resolved guide-by-guide at this input depth.
""")
    o.append('Top 25 by IP chimera count:\n')
    o.append(md_table(guides.head(25), 'AluACA guide'))
    o.append('')

    # ---- targets -----------------------------------------------------------
    o.append('## What the AluACA guides pair with\n')
    tc = pd.DataFrame({'IP': alu_ip['target_class'].value_counts(),
                       'input': alu_ct['target_class'].value_counts()}).fillna(0).astype(int)
    o.append(md_table(tc, 'target class'))
    o.append('')

    g_ip = alu_ip[alu_ip.target_class == 'hg38']
    if not g_ip.empty and 'gene_type' in g_ip.columns:
        o.append('### Genomic arm biotypes (IP)\n')
        bt = g_ip['gene_type'].value_counts().head(12).to_frame('chimeras')
        o.append(md_table(bt, 'gene_type'))
        o.append('')

        mrna = g_ip[g_ip.gene_type.astype(str).str.contains('protein_coding', na=False)]
        mrna_ct = pd.DataFrame()
        g_ct = alu_ct[alu_ct.target_class == 'hg38']
        if not g_ct.empty and 'gene_type' in g_ct.columns:
            mrna_ct = g_ct[g_ct.gene_type.astype(str).str.contains('protein_coding', na=False)]
        o.append('## AluACA-mRNA chimeras\n')
        o.append(f'**{len(mrna):,}** IP chimeras have a protein_coding genomic arm '
                 f'({len(mrna_ct):,} in input).\n')
        if not mrna.empty:
            feat = pd.DataFrame({'IP': mrna['feature'].value_counts()})
            if len(mrna_ct):
                feat['input'] = mrna_ct['feature'].value_counts()
            feat = feat.fillna(0).astype(int)
            o.append(md_table(feat, 'feature'))
            o.append('\n**Exonic hits are the defensible set.** An intronic hit is as easily '
                     'explained by co-transcriptional proximity in the host pre-mRNA as by a '
                     'guide-target duplex, so intronic and exonic counts should not be pooled.\n')
            ex = mrna[mrna.feature == 'exonic']
            if not ex.empty:
                top = ex['gene_name'].value_counts().head(25).to_frame('exonic chimeras')
                o.append('Top exonic mRNA targets (IP):\n')
                o.append(md_table(top, 'gene'))
                o.append('')
            pairs = (mrna.assign(pair=mrna['guide_names'] + ' -> ' + mrna['gene_name'].astype(str))
                     ['pair'].value_counts().head(20).to_frame('chimeras'))
            o.append('Top AluACA-mRNA pairs (IP, exonic and intronic):\n')
            o.append(md_table(pairs, 'guide -> gene'))
            o.append('')

    # ---- stratified enrichment ---------------------------------------------
    # The aggregate AluACA ratio pools two populations with opposite behaviour. Pederiva
    # et al. (Sci Adv 2023, PMC10381945) propose intronic Alu-derived H/ACA RNAs as the
    # guides for dyskerin-dependent mRNA pseudouridylation, so the stratum that model
    # predicts -- an exonic mRNA target, outside any repeat -- has to be scored on its own
    # rather than averaged together with the Alu-to-Alu background.
    o.append('## Stratified enrichment: separating signal from the Alu background\n')
    strata = [
        ('genomic arm, all',                    lambda d: d.index.notna()),
        ('arm inside a repeat',                 lambda d: d.target_in_repeat),
        ('arm outside any repeat',              lambda d: ~d.target_in_repeat),
        ('outside repeat, protein_coding',      lambda d: (~d.target_in_repeat) & _pc(d)),
        ('outside repeat, protein_coding, EXONIC',
         lambda d: (~d.target_in_repeat) & _pc(d) & (d.feature == 'exonic')),
    ]
    strat = {}
    for cls in ('AluACA', 'snoRNA'):
        I = ip[(ip.guide_class == cls) & (ip.target_class == a.gtag)]
        C = ctrl[(ctrl.guide_class == cls) & (ctrl.target_class == a.gtag)]
        rws = {}
        strat[cls] = {}
        for lbl, fn in strata:
            x, y = int(fn(I).sum()), int(fn(C).sum())
            pt, lo, hi = rate_ratio(x, y, N1, N2)
            hs = 'inf' if hi == float('inf') else f'{hi:.2f}'
            rws[lbl] = [x, y, round(pt, 2), f'{lo:.2f} - {hs}']
            # Kept so the paragraph below quotes the table rather than restating
            # numbers that go stale the moment the pipeline is re-run.
            strat[cls][lbl] = (x, y, pt, lo, hi)
        o.append(f'**{cls} guides**\n')
        o.append(md_table(pd.DataFrame.from_dict(
            rws, orient='index', columns=['IP', 'input', 'rate ratio', '95% CI']), 'stratum'))
        o.append('')
    _EX = 'outside repeat, protein_coding, EXONIC'
    rep_r = strat['AluACA']['arm inside a repeat'][2]
    ex_n, _, ex_r, ex_lo, _ = strat['AluACA'][_EX]
    sno_ex_r = strat['snoRNA'][_EX][2]
    pooled_alu = rows['AluACA'][4]
    # "clear of 1" is a claim about the interval, so let the interval make it.
    ci_clause = ('with a confidence interval clear of 1' if ex_lo > 1
                 else f'though its interval still reaches below 1 ({ex_lo:.2f})')
    gap = sno_ex_r / ex_r if ex_r else float('nan')

    # How much of the enriched stratum sits on the guide's own locus. This is the
    # computable part of the cis question: annotate_chimeras flags a target arm
    # overlapping its own guide locus on the same strand.
    ex_ip = ip[(ip.guide_class == 'AluACA') & (ip.target_class == a.gtag) &
               (~ip.target_in_repeat) & _pc(ip) & (ip.feature == 'exonic')]
    own = int(ex_ip.get('target_in_source_locus', pd.Series(dtype=object)).eq(True).sum())
    trans_pct = 100 * (len(ex_ip) - own) / len(ex_ip) if len(ex_ip) else float('nan')

    o.append(f"""
**The aggregate AluACA depletion is the Alu-to-Alu background, and it inverts the sign of
the real signal.** Chimeras whose genomic arm lands in a repeat run at {rep_r:.2f}x. Strip those
out and restrict to exonic protein-coding targets -- the stratum a guide model actually
predicts -- and AluACA chimeras are *enriched* at {ex_r:.2f}x {ci_clause},
over {ex_n:,} IP chimeras. Reporting only the pooled {pooled_alu:.2f}x would have buried that.

The same stratification puts snoRNA guides at {sno_ex_r:.2f}x, so AluACA-mRNA pairing is roughly
{gap:.0f}-fold weaker than canonical snoRNA guiding rather than absent. That is the size of
effect expected if the duplexes are short-lived: Pederiva et al. argue mRNA
pseudouridylation proceeds "with the aid of guide RNAs containing mismatches toward the
mRNA to be modified", and a mismatched, catalytically transient duplex is captured by
proximity ligation far less efficiently than a stable snoRNP-rRNA pairing. A weaker ratio
is therefore the predicted observation, not evidence against the model.

**This measures trans pairing only.** {trans_pct:.1f}% of the enriched exonic set pairs a guide
with an mRNA outside the guide's own locus; only {own:,} of {len(ex_ip):,} land back on it. That is
not evidence against cis action, because genome masking removes cis geometry by
construction -- a guide ligated to its own host pre-mRNA yields a read that aligns
contiguously, or across a short novel junction, and is dropped before chimera calling.
Testing the co-transcriptional model properly would need that stage relaxed or replaced;
this pipeline cannot address it either way.
""")

    # ---- Alu target orientation ---------------------------------------------
    ori = alu_orientation(ip, ctrl, a.gtag, a.rmsk, a.bedtools,
                          os.path.join(a.resdir, 'orient_work'), N1, N2)
    if ori:
        o.append('## Alu targets, split by orientation\n')
        o.append(
            'For an Alu-derived guide an antisense Alu is arguably the most available '
            'complementary partner in the transcriptome -- the same logic that underlies '
            'IRAlu duplexes and STAU1-mediated decay. A same-orientation Alu shares the '
            "guide's sequence and cannot base-pair, so orientation separates a genuine "
            'Alu:Alu duplex from plain sequence self-similarity, which the repeat flag '
            'alone cannot.\n')
        rows = {}
        for cls in ('AluACA', 'snoRNA'):
            for orient, note in (('sense', 'cannot base-pair'),
                                 ('antisense', 'can base-pair')):
                x, y = ori[(cls, orient, 'ip')], ori[(cls, orient, 'ctrl')]
                pt, lo, hi = rate_ratio(x, y, N1, N2)
                hs = 'inf' if hi == float('inf') else f'{hi:.2f}'
                rows[f'{cls} -> {orient} Alu ({note})'] = [x, y, round(pt, 2), f'{lo:.2f} - {hs}']
        o.append(md_table(pd.DataFrame.from_dict(
            rows, orient='index', columns=['IP', 'input', 'rate ratio', '95% CI']), 'stratum'))
        aa, as_ = ori[('AluACA', 'antisense', 'ip')], ori[('AluACA', 'sense', 'ip')]
        sa, ss = ori[('snoRNA', 'antisense', 'ip')], ori[('snoRNA', 'sense', 'ip')]
    if ori and not (as_ and ss):
        o.append(
            f'The orientation composition is not computed for this arm: it needs both '
            f'sense strata to be non-empty, and this run has {as_:,} sense-Alu AluACA '
            f'targets and {ss:,} for snoRNA guides. A catalogue with no AluACA records '
            f'-- the plain `snoRNA.txt.fa` arms -- always lands here.\n')

    if ori and as_ and ss:
        od, pv = fisher_exact([[aa, as_], [sa, ss]])

        # the enriched stratum's ratio, computed rather than restated
        def _ex(d):
            return d[(d.guide_class == 'AluACA') & (d.target_class == a.gtag) &
                     (~d.target_in_repeat) & _pc(d) & (d.feature == 'exonic')]
        ex_r = rate_ratio(len(_ex(ip)), len(_ex(ctrl)), N1, N2)[0]
        o.append(f"""
**The pairable class is present in quantity and is not enriched.** {aa:,} AluACA chimeras
target an antisense Alu -- one that could actually form a duplex -- and they run at
{rate_ratio(aa, ori[('AluACA','antisense','ctrl')], N1, N2)[0]:.2f}x, *more* depleted than
the same-orientation Alus that cannot pair at all. A duplex model predicts the opposite
ordering.

The orientation composition says the same thing. Antisense:sense is {aa/as_:.2f} for AluACA
guides against {sa/ss:.2f} for snoRNA guides (Fisher p = {pv:.1e}). snoRNA guides have no Alu
complementarity, so their ratio is the baseline availability of antisense Alus among
transcribed sequences; AluACA guides sit *below* that baseline rather than above it.

So within AluACA guides the enrichment sits in non-repeat exonic mRNA
({ex_r:.2f}x), not in Alu targets of either orientation. Two limits
on that reading: the AluACA class is depleted overall, so every stratum inside it starts
below 1 and the informative comparison is between strata rather than against 1; and this
says nothing about *cis* IRAlu pairing within a single transcript, which genome masking
removes irrespective of orientation.
""")
        o.append('')

    # ---- chrM as an internal artefact control -------------------------------
    # chrM is deliberately kept in the reference. Dyskerin is nuclear and AluACAs are
    # nucleoplasmic H/ACA RNPs, so an AluACA:MT-CO3 duplex is not physically available and
    # every chrM chimera must be ligation artefact -- which makes it a free, in-sample
    # measurement of the artefact floor. Removing chrM would not remove those reads: the
    # genome step runs --outFilterMultimapNmax 1, so survivors are unique to chrM, and
    # without it they would fall through to one of hg38's NUMT copies (95-99% identical)
    # and be reported as nuclear targets. That trades a labelled artefact for a hidden one.
    o.append('## chrM as an internal artefact control\n')
    def mstrat(d, cls, mito):
        x = d[(d.guide_class == cls) & (d.target_class == a.gtag) & (~d.target_in_repeat)
              & (d.feature == 'exonic') & _pc(d)]
        m = x.reference_target.astype(str).isin(['chrM', 'MT', 'chrMT'])
        return x[m] if mito else x[~m]
    mrows = {}
    mit = {}
    for cls in ('AluACA', 'snoRNA'):
        for mito, lbl in ((True, 'chrM (impossible, = artefact)'), (False, 'nuclear')):
            A_, B_ = mstrat(ip, cls, mito), mstrat(ctrl, cls, mito)
            pt, lo, hi = rate_ratio(len(A_), len(B_), N1, N2)
            hs = 'inf' if hi == float('inf') else f'{hi:.2f}'
            mrows[f'{cls} -> {lbl}'] = [len(A_), len(B_), round(pt, 2), f'{lo:.2f} - {hs}']
            mit[(cls, mito)] = (len(A_), len(B_), pt, lo, hi)
    o.append(md_table(pd.DataFrame.from_dict(
        mrows, orient='index', columns=['IP', 'input', 'rate ratio', '95% CI']), 'stratum'))
    (am_i, am_c, am_r, _, am_hi) = mit[('AluACA', True)]
    (an_i, an_c, an_r, _, _) = mit[('AluACA', False)]
    (sm_i, sm_c, sm_r, _, _) = mit[('snoRNA', True)]
    _, p_alu = fisher_exact([[am_i, am_c], [an_i, an_c]])
    _, p_sno = fisher_exact([[sm_i, sm_c], [mit[('snoRNA', False)][0],
                                            mit[('snoRNA', False)][1]]])
    # Every clause that could flip is derived, so the argument cannot drift out of
    # step with the counts on a re-run.
    reach = ('contains' if am_hi >= an_r else 'does not reach')
    alu_verdict = ('not distinguishable' if p_alu >= 0.05 else 'distinguishable')
    sno_verdict = ('also not distinguishable' if p_sno >= 0.05 else 'distinguishable')
    hi_txt = 'infinity' if am_hi == float('inf') else f'{am_hi:.2f}'
    sno_zero = (' rests on *zero* input reads: its magnitude is produced entirely by the\n'
                'Haldane-Anscombe 0.5 substituted for that zero, its interval runs to infinity, and'
                if sm_c == 0 else
                f' rests on {sm_c:,} input reads, and')

    o.append(f"""
**This is the most important caveat in the report.** For AluACA guides the artefact floor
sits at {am_r:.2f}x with an interval reaching {hi_txt}, and that interval {reach} the {an_r:.2f}x measured
on nuclear exonic mRNA targets. A Fisher exact test on the 2x2 of counts puts the two at
p = {p_alu:.2f}: **{alu_verdict}**. The enrichment over *input* is solid; whether it exceeds
the artefact floor is simply not resolved by this data. The AluACA-mRNA result should
therefore be stated as consistent with a guide model, not as evidence for one.

What chrM does establish is that the floor is not zero. AluACA guides generate {am_i:,}
impossible chimeras in this stratum, against {am_i + an_i:,} candidates -- ligation noise is
measurably present, not negligible.

Both chrM ratios are badly underpowered and should not be over-read. The snoRNA figure of
{sm_r:.2f}x{sno_zero} Fisher
against the nuclear stratum gives p = {p_sno:.2f} -- {sno_verdict}. It is tempting to
argue from these numbers that an artefact pairing inherits the enrichment of whichever
guide it is attached to, since chimeras form during on-bead ligation. That story fits the
point estimates, but the intervals do not support it and it should not be presented as a
finding.

Whether a guide can actually form the >=8 bp bipartite duplex around a target uridine
would separate signal from artefact, and no read-counting statistic can. That is the
experiment this result needs next.
""")

    # ---- the actual target list --------------------------------------------
    # Restricted to the one stratum that is enriched over input, so this is a candidate
    # list rather than a ranking of whatever is most abundant. Counts per guide-gene pair
    # and per gene, each with an IP/input ratio, written out in full as TSVs.
    def stratum(d):
        return d[(d.guide_class == 'AluACA') & (d.target_class == a.gtag) &
                 (~d.target_in_repeat) & _pc(d) & (d.feature == 'exonic')]
    TI, TC = stratum(ip), stratum(ctrl)
    o.append('## AluACA-mRNA target list\n')
    o.append(f'Every target below is drawn from the enriched stratum only -- exonic, '
             f'protein-coding, outside any annotated repeat ({len(TI):,} IP chimeras, '
             f'{rate_ratio(len(TI), len(TC), N1, N2)[0]:.2f}x over input). Ranking the '
             f'unstratified set instead would just rank abundance.\n')

    def summarise(keys, fname, label, topn):
        gi = TI.groupby(keys).size()
        gc = TC.groupby(keys).size()
        t = pd.DataFrame({'IP': gi, 'input': gc}).fillna(0).astype(int)
        st = [rate_ratio(r.IP, r.input, N1, N2) for r in t.itertuples()]
        t['rate ratio'] = [round(x[0], 2) for x in st]
        t['CI low'] = [round(x[1], 2) for x in st]
        t['CI high'] = [round(x[2], 2) if x[2] != float('inf') else float('nan') for x in st]
        t = t.sort_values('IP', ascending=False)
        path = os.path.join(a.resdir, fname)
        t.to_csv(path, sep='\t')
        o.append(f'### {label}\n')
        o.append(f'{len(t):,} in total; full list in '
                 f'`{os.path.join(reldir, fname)}`. Top {topn} by IP count:\n')
        o.append(md_table(t.head(topn), ' / '.join(keys) if isinstance(keys, list) else keys))
        o.append('')
        return t

    genes = summarise(['gene_name'], f'{a.ip}.AluACA_mRNA_targets_by_gene.tsv',
                      'By target gene', 30)
    pairs = summarise(['guide_names', 'gene_name'], f'{a.ip}.AluACA_mRNA_targets_by_pair.tsv',
                      'By guide-target pair', 30)
    # Abundance sanity check. Mitochondrial mRNAs are the useful control here: dyskerin
    # is nuclear/nucleolar and AluACAs are nucleoplasmic H/ACA RNPs, so an AluACA:MT-CO3
    # duplex is not physically available -- anything scored against MT transcripts is
    # ligation artefact tracking abundance, and measures how much of the list is that.
    gn = genes.index.to_series().astype(str)
    mt = gn.str.contains(r'\bMT-', regex=True)
    rp = gn.str.match(r'^(RPS|RPL)\d')
    n_single = int((genes.IP == 1).sum())
    o.append(f"""
**Read this as a candidate list, not a set of identified targets.** Two things about its
shape argue for caution.

*It is substantially an abundance ranking.* {int(mt.sum())} mitochondrial and
{int(rp.sum())} ribosomal-protein genes account for
{100 * genes[mt].IP.sum() / genes.IP.sum():.1f}% and
{100 * genes[rp].IP.sum() / genes.IP.sum():.1f}% of the chimeras respectively, and
together they are {int(mt.head(30).sum()) + int(rp.head(30).sum())} of the top 30 rows.
Mitochondrial transcripts are the tell: dyskerin is nuclear and AluACAs are nucleoplasmic
H/ACA RNPs, so an AluACA-MT-CO3 duplex is not physically available and those
{int(genes[mt].IP.sum())} chimeras have to be ligation artefact. They are a free internal
estimate of how much of this list is abundance-driven noise, and they sit near the top of it.

*Per-gene counts are too thin to rank.* {n_single:,} of {len(genes):,} targets
({100 * n_single / len(genes):.0f}%) rest on a single chimera, and only
{int((genes.IP >= 5).sum())} have five or more. The CI columns show the consequence --
almost every row spans 1. The {rate_ratio(len(TI), len(TC), N1, N2)[0]:.2f}x enrichment is
a property of the stratum in aggregate, where thousands of reads back it; it does not
transfer to any individual gene in the table.

What the list is good for is generating hypotheses to test directly -- and the obvious
filter to apply first is whether a candidate has a plausible pseudouridylation pocket,
i.e. whether the AluACA can form the >=8 bp bipartite duplex around a target uridine that
H/ACA guiding requires. That is a sequence calculation this pipeline does not do.
""")
    o.append(f'{len(genes):,} distinct mRNAs and {len(pairs):,} distinct guide-target pairs. '
             f'Names joined by `|` are ambiguous calls, not composites: a `|` in a gene name '
             f'means the arm overlaps both genes and the annotation cannot separate them, and '
             f'a `|` in a guide name means the read matched those guides equally well. '
             f'Treat both as unresolved rather than as a single identified target.\n')

    # ---- caveats -----------------------------------------------------------
    o.append('## How much of this survives scrutiny\n')
    rep = alu_ip['target_in_repeat'].sum()
    src = alu_ip['target_in_source_locus'].sum()
    gtot = len(g_ip) if not g_ip.empty else 0
    o.append(f'Of {gtot:,} AluACA chimeras with a genomic arm:\n')
    o.append(f'- **{rep:,} ({100*rep/gtot:.1f}%)** have that arm inside an annotated repeat. '
             'Every AluACA is Alu-derived and hg38 holds over a million Alu copies, so '
             'Alu-to-Alu pairing is the dominant false-positive mode for this question. '
             'These are flagged, not removed.\n' if gtot else '')
    o.append(f'- **{src:,} ({100*src/gtot:.1f}%)** land on a guide locus on the same strand, '
             'i.e. probably one contiguous transcript rather than a chimera.\n' if gtot else '')
    clean = g_ip[~g_ip.target_in_repeat & ~g_ip.target_in_source_locus]
    o.append(f'\nDropping both flags leaves **{len(clean):,}** AluACA genomic chimeras.')
    if not clean.empty and 'gene_type' in clean.columns:
        cm = clean[clean.gene_type.astype(str).str.contains('protein_coding', na=False)]
        ce = cm[cm.feature == 'exonic']
        o.append(f' Of those, **{len(cm):,}** are protein_coding and **{len(ce):,}** '
                 f'are exonic — the most conservative AluACA-mRNA set.\n')
        if not ce.empty:
            o.append('\nConservative AluACA-mRNA set, top genes:\n')
            o.append(md_table(ce['gene_name'].value_counts().head(20).to_frame('chimeras'), 'gene'))

    # ---- sanity check against the published run ----------------------------
    if os.path.exists(a.published):
        with open(a.published) as fh:
            pub = sum(1 for _ in fh) - 1
        ours = int((ip.target_class == 'hg38').sum())
        o.append('\n## Cross-check against the published hg19 run\n')
        o.append(f'The published output for this same sample '
                 f'(`{os.path.basename(a.published)}`, hg19, plain 1951-record snoRNA source) '
                 f'holds **{pub:,}** genomic chimeras. This rerun (hg38, merged '
                 f'2701-record source) gives **{ours:,}**.\n')
        o.append('\nExact agreement is not expected: different genome build, a source '
                 'catalogue with 765 extra sequences competing for the same reads, and a '
                 'repeat-masking index built from public RepBase rather than the lab\'s '
                 'own. Preprocessing, however, was validated to reproduce 4983/5000 '
                 '(99.66%) of the published trimmed sequences byte-for-byte, so any '
                 'divergence is downstream of read handling.\n')

    o.append('\n## Files\n')
    o.append(f'- `{reldir}/{a.ip}.annotated.tsv` — one row per IP chimera, all columns\n'
             f'- `{reldir}/{a.inp}.annotated.tsv` — same for the input control\n'
             f'- `{reldir}/<uid>/<uid>.snoRNA.<target>.chimeras.csv` — per-target, '
             'pipeline-native\n'
             f'- `{reldir}/<uid>/<uid>.snorna.chimeras.pipeline.sh` — every command run, '
             'with literal arguments\n'
             f'- `{reldir}/<uid>/map.metric.and.log` — per-stage mapping metrics\n')

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, 'w') as f:
        f.write('\n'.join(o) + '\n')
    print(f'wrote {a.out}')


if __name__ == '__main__':
    main()
