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
from scipy.stats import beta, chi2, fisher_exact

# paths.py lives one level up, shared with the analysis scripts; the repo is a
# collection of scripts rather than an installed package, so put src/ on the path.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import find_input, find_tool, proj
from annotate_chimeras import MIN_TARGET_ARM

CONF = 0.95
SRC_NOTE = 'Song et al. 2025, Genome Biology, doi:10.1186/s13059-025-03508-7'
MANAKOV = 'Manakov et al. 2022, bioRxiv, doi:10.1101/2022.02.13.480296'
# A guide class, stratum or guide is only said to carry real chimeras when the upper
# bound of its false-positive share is below this. 1.0 = "input cannot explain all of it".
SIGNAL = 1.0


def share_ci(x, y, n1, n2, conf=CONF):
    """Estimated false-positive share of the IP calls, (y/n2) / (x/n1), with an exact CI.

    The input control is assumed not to go through the chimeric ligation (the protocol does
    not say either way), so every chimera called in it is treated as a pipeline or
    library-prep artefact. Its per-read rate, set against the IP's,
    estimates what share of the IP calls those same artefacts account for -- on the
    assumption that artefacts arise at the same rate per trimmed read in both libraries.

    x ~ Pois(l1*n1) and y ~ Pois(l2*n2). Conditional on n = x+y, y ~ Binomial(n, p) with
    share = l2/l1 = (p/(1-p)) * (n1/n2), so a Clopper-Pearson interval on p maps
    monotonically onto the share. No pseudocount: y = 0 is a point estimate of 0 with a
    finite upper bound, and x = 0 has no share at all.
    """
    if x == 0:
        return float('nan'), float('nan'), float('nan')
    n = x + y
    a = (1 - conf) / 2
    p_lo = 0.0 if y == 0 else beta.ppf(a, y, n - y + 1)
    p_hi = beta.ppf(1 - a, y + 1, n - y)
    f = lambda q: (q / (1 - q)) * (n1 / n2)
    return (y / n2) / (x / n1), f(p_lo), f(p_hi)


def poisson_ci(k, conf=CONF):
    """Exact (Garwood) interval on a Poisson count; k = 0 still has an upper bound."""
    a = (1 - conf) / 2
    lo = 0.0 if k == 0 else chi2.ppf(a, 2 * k) / 2
    return lo, chi2.ppf(1 - a, 2 * k + 2) / 2


def verdict(hi):
    """What the upper bound of the false-positive share supports -- and no more.

    The input cannot see artefactual chimeras that form after lysis, so a call the input
    does not explain is IP-specific: necessary for real pairing, not sufficient. The
    labels say that rather than "real"."""
    if hi != hi:  # NaN: no IP calls
        return 'no IP calls'
    if hi < 0.5:
        return 'mostly IP-specific'
    if hi < SIGNAL:
        return 'partly IP-specific'
    return 'not distinguishable from input artefacts'


def _pct(v):
    return 'inf' if v == float('inf') else f'{100 * v:.1f}%'


def fp_cells(x, y, n1, n2):
    """The standard columns for one row: counts, input rate, share, estimated real calls."""
    sh, lo, hi = share_ci(x, y, n1, n2)
    r_lo, r_hi = poisson_ci(y)
    if x:
        real = f'{max(0, x * (1 - sh)):,.0f} ({max(0, x * (1 - hi)):,.0f} - {max(0, x * (1 - lo)):,.0f})'
        share = f'{_pct(sh)} ({_pct(lo)} - {_pct(hi)})'
    else:
        real = share = ''
    return [x, y, round(1e6 * x / n1, 1),
            f'{1e6 * y / n2:.1f} ({1e6 * r_lo / n2:.1f} - {1e6 * r_hi / n2:.1f})',
            share, real, verdict(hi)]


FP_COLS = ['IP', 'input', 'IP per M', 'input per M (95% CI)',
           'false-positive share (95% CI)', 'est. real IP chimeras', 'verdict']

FLAGS = ('guide_low_complexity', 'target_arm_short', 'target_arm_aligned',
         'read_contiguous', 'guide_near_target', 'genome_contiguity')


def _flag(d, c):
    return d[c].astype(str).str.lower().eq('true')


def usable(d, gtag):
    """Calls that survive the artefact flags from annotate_chimeras.py.

    For every target: the guide arm must not be a simple repeat (DUST), the target arm
    must not be too short to place (`target_arm_short`), and the whole read must not align
    contiguously anywhere in the genome, nor map to too many loci to check. A genomic call must also re-align at its reported locus (otherwise it
    cannot be checked) and must not have its guide arm within 2 kb of the target."""
    genomic = d.target_class == gtag
    ok = ~_flag(d, 'target_arm_aligned') | _flag(d, 'guide_near_target')
    gw = d.genome_contiguity.isin(['contiguous', 'too many loci'])
    return (~_flag(d, 'guide_low_complexity') & ~_flag(d, 'target_arm_short')
            & ~gw & ~(genomic & ok))


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
    missing = [c for c in FLAGS if c not in df.columns]
    if not missing and df['genome_contiguity'].isna().all():
        missing = ['genome_contiguity (column present but empty: no STAR index was given)']
    if missing:
        sys.exit(f'{path} has no artefact flags ({", ".join(missing)}). Add them with\n'
                 f'  python3 src/chimeric/annotate_chimeras.py --annotated {path} --out {path} '
                 f'--gtag <build> --genome-index <STAR index>\n'
                 f'The genome-wide check loads the STAR index, so run it where that fits in RAM.')
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


def crosstab(df, label, gtag='hg38'):
    ct = pd.crosstab(df['guide_class'], df['target_class'])
    for c in ('rRNA', 'snRNA', 'tRNA', gtag):
        if c not in ct.columns:
            ct[c] = 0
    ct = ct[[gtag, 'rRNA', 'snRNA', 'tRNA']]
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
    o.append(f'# AluACA-guided chimeric reads in DKC1 chimeric eCLIP ({a.gtag})\n')
    o.append('Generated by `src/chimeric/make_report.py` from the pipeline output; '
             're-run it after any rerun rather than editing numbers here. '
             'The procedure, the reference build and every deviation from upstream '
             'are documented in [`src/chimeric/README.md`](../../src/chimeric/README.md).\n')

    o.append('## Samples\n')
    samp = pd.DataFrame({
        'role': ['IP', 'input (assumed unligated: false-positive control)'],
        'GSM': ['GSM8521923', 'GSM8521922'],
        'raw reads': [n_ip.get('raw', 0), n_ct.get('raw', 0)],
        'after trimming': [n_ip.get('trimmed', 0), n_ct.get('trimmed', 0)],
        'after repeat+genome masking': [n_ip.get('after_genome_mask', 0),
                                        n_ct.get('after_genome_mask', 0)],
    }, index=[a.ip, a.inp])
    o.append(md_table(samp, 'run'))
    o.append('\nMasking keeps only reads that fail to align end-to-end to both the RepBase '
             f'human consensus set and {a.gtag} — the pipeline\'s definition of a candidate '
             'chimera. It is not a complete filter: a contiguous read with a few non-genomic '
             'bases at its ends (a poly(A) tail, an adapter remnant), or one ~95% identical '
             'to several repeat copies, also fails end-to-end alignment and goes on to chimera '
             'calling. That is what the artefact flags below catch.\n')

    o.append('## Chimeras by guide class and target\n')
    o.append(crosstab(ip, f'{a.ip} — IP', a.gtag))
    o.append('')
    o.append(crosstab(ctrl, f'{a.inp} — input', a.gtag))
    o.append('\n`guide_class` is decided by exact membership of the 765 AluACA union names, '
             'not by identifier pattern: 92 names in `snoRNA.txt.fa` also begin `id3` and a '
             'prefix test miscounts them as AluACA. `ambiguous` means a read\'s equal-scoring '
             'guides span both catalogues.\n')

    # ---- what the input measures, and the flags ------------------------------
    N1 = n_ip.get('trimmed') or 1
    N2 = n_ct.get('trimmed') or 1
    ip['usable'] = usable(ip, a.gtag)
    ctrl['usable'] = usable(ctrl, a.gtag)
    U_ip, U_ct = ip[ip.usable], ctrl[ctrl.usable]
    classes = [c for c in ('snoRNA', 'AluACA', 'ambiguous') if (ip.guide_class == c).any()]
    have_alu = 'AluACA' in classes

    o.append('## What the input control measures\n')
    rr_ip = int(((U_ip.guide_class == 'snoRNA') & (U_ip.target_class == 'rRNA')).sum())
    rr_ct = int(((U_ct.guide_class == 'snoRNA') & (U_ct.target_class == 'rRNA')).sum())
    rr_ipm, rr_ctm = 1e6 * rr_ip / N1, 1e6 * rr_ct / N2
    rr_fold = (f'a {rr_ipm / rr_ctm:,.0f}-fold difference' if rr_ct
               else 'with none at all in the input')
    o.append(f"""The chimeric ligation is performed on the beads after immunoprecipitation, and the input
control is 2% of the sample saved before the IP ({SRC_NOTE}, methods). Neither that paper
nor the GEO protocol says whether the input then goes through the chimeric ligation.
**This report assumes it does not**, and treats every chimera called in the input as a
false positive from the pipeline or from library preparation.

The data are consistent with that assumption but do not prove it. snoRNA-guided rRNA
chimeras, the canonical DKC1 pairing, occur at {rr_ipm:,.1f} usable calls per million trimmed
reads in the IP against {rr_ctm:,.1f} in the input, {rr_fold}. Ligation of a dilute lysate
that simply worked badly would also give few chimeras. If the assumption is wrong, the
input contains some real chimeras, and every false-positive share below is an
overestimate -- the direction that makes a signal harder, not easier, to find.

Under the assumption, the input does not give a background *rate of chimera formation*
to divide by; it gives a *false-positive rate* for the calling procedure.

This report therefore does two things. It removes the artefacts the input calls turned
out to be (the flags from `annotate_chimeras.py`, below), and it uses what is left in the
input to estimate the **false-positive share** of the IP calls: the input's calls per
million trimmed reads divided by the IP's. A share whose 95% upper bound is below 100%
means the input cannot account for all of the IP calls: some are **IP-specific**. That is
necessary for real pairing but not sufficient, because artefactual chimeras that form
after lysis are IP-specific too (see chrM). A share whose upper bound is at or above 100% means the
data cannot distinguish that class from the artefacts the input does measure.

Arithmetically the share is the inverse of an IP/input rate ratio. What changes is the
reading, that it is computed after the artefact flags, and that a zero input count is
reported as an upper bound instead of being replaced by 0.5.

Two assumptions limit it. Artefacts must arise at the same rate per trimmed read in both
libraries, and the libraries did not trim alike ({100 * n_ip.get('trimmed', 0) / max(n_ip.get('raw', 1), 1):.0f}% of IP reads survived trimming,
{100 * n_ct.get('trimmed', 0) / max(n_ct.get('raw', 1), 1):.0f}% of input reads), so the share is an estimate, not a measurement. And, under the
assumption above, the input cannot measure chimeras that form after lysis between RNAs
that were not paired in the cell. Those are documented for chimeric eCLIP: mixing human
and rat lysates before an AGO2 IP gave 8.6% of human miRNA chimeras with rat targets,
against a 1.2% baseline, and diluting the beads lowered the rate without reaching
significance, so crowding on the beads and complexes meeting new RNA after lysis were not
separated ({MANAKOV}). That rate was measured for AGO2 and may differ for DKC1. The part
of it this data can see is estimated from chrM, below.
""")

    o.append('### Artefact flags\n')
    o.append('Each call is attributed to the first flag that removes it, so the rows sum. '
             'A simple-repeat guide arm matches too many sequences to say which RNA it came '
             'from. "Contiguous in the genome" means one alignment of the whole read, at any '
             'locus, covers both arms: a single transcript split in two. The last two columns '
             'apply to genomic targets only.\n')
    frows = {}
    for lib, d in (('IP', ip), ('input', ctrl)):
        for cls in classes:
            x = d[d.guide_class == cls]
            lowc = _flag(x, 'guide_low_complexity')
            gen = x.target_class == a.gtag
            short = ~lowc & _flag(x, 'target_arm_short')
            cont = ~lowc & ~short & x.genome_contiguity.eq('contiguous')
            many = ~lowc & ~short & ~cont & x.genome_contiguity.eq('too many loci')
            rest = ~lowc & ~short & ~cont & ~many
            unal = gen & rest & ~_flag(x, 'target_arm_aligned')
            near = gen & rest & ~unal & _flag(x, 'guide_near_target')
            frows[f'{lib} {cls}'] = [len(x), int(lowc.sum()), int(short.sum()), int(cont.sum()),
                                     int(many.sum()), int(unal.sum()), int(near.sum()),
                                     int(x.usable.sum()),
                                     f'{100 * x.usable.mean():.1f}%' if len(x) else '']
    o.append(md_table(pd.DataFrame.from_dict(frows, orient='index', columns=[
        'called', 'simple-repeat guide (DUST)', 'short target arm', 'contiguous in the genome',
        'too many loci', 'genomic, not re-aligned', 'genomic, guide within 2 kb', 'usable',
        'usable share']), 'library / guide class'))
    o.append('')

    # ---- how much the length threshold matters -------------------------------
    # The threshold is a judgement call, so print the whole series rather than only the
    # chosen value: a reader can see whether a conclusion depends on where it sits.
    o.append('### Sensitivity to the target-arm length threshold\n')
    base_ip = ip[ip.usable | _flag(ip, 'target_arm_short')]
    base_ct = ctrl[ctrl.usable | _flag(ctrl, 'target_arm_short')]
    lens = [16, 20, 25, 30, 35, 40]
    srows = {}
    for cls in [c for c in ('snoRNA', 'AluACA') if c in classes]:
        for m in lens:
            I = base_ip[(base_ip.guide_class == cls) & (base_ip.target_class == a.gtag)
                        & (base_ip.map_to_target_length >= m)]
            C = base_ct[(base_ct.guide_class == cls) & (base_ct.target_class == a.gtag)
                        & (base_ct.map_to_target_length >= m)]
            sh, lo, hi = share_ci(len(I), len(C), N1, N2)
            srows[f'{cls}, arm >= {m} nt'] = [
                len(I), len(C), round(1e6 * len(I) / N1, 1), round(1e6 * len(C) / N2, 1),
                f'{_pct(sh)} ({_pct(lo)} - {_pct(hi)})' if len(I) else '', verdict(hi)]
    o.append(md_table(pd.DataFrame.from_dict(srows, orient='index', columns=[
        'IP', 'input', 'IP per M', 'input per M', 'false-positive share (95% CI)', 'verdict']),
        'guide class, minimum target arm'))
    o.append(f"""
Every other flag is applied here; only the length cut varies. The report's tables use
{MIN_TARGET_ARM} nt. The pipeline itself accepts 16 nt, where a unique placement in a 3 Gb genome
is largely chance -- and the input's calls sit at that end, so raising the cut removes
input calls much faster than IP calls. A conclusion that appears only below the chosen
threshold is a conclusion about short arms, not about pairing.
""")

    # ---- false-positive share by guide class ---------------------------------
    o.append('## False-positive estimate by guide class\n')
    groups = [(a.gtag, 'genome', lambda d: d.target_class == a.gtag),
              ('rRNA', 'rRNA', lambda d: d.target_class == 'rRNA'),
              ('RNA', 'snRNA + tRNA', lambda d: d.target_class.isin(['snRNA', 'tRNA']))]
    rows, res = {}, {}
    for cls in classes:
        for key, lbl, fn in groups:
            I, C = U_ip[(U_ip.guide_class == cls) & fn(U_ip)], U_ct[(U_ct.guide_class == cls) & fn(U_ct)]
            if len(I) == 0 and len(C) == 0:
                continue
            rows[f'{cls} -> {lbl}'] = fp_cells(len(I), len(C), N1, N2)
            res[(cls, key)] = (len(I), len(C)) + share_ci(len(I), len(C), N1, N2)
    o.append(md_table(pd.DataFrame.from_dict(rows, orient='index', columns=FP_COLS),
                      'guide -> target'))
    # A share whose lower bound exceeds 100% cannot arise under the equal-artefact-rate
    # assumption, so it is evidence against that assumption for the class.
    alu_uncal = ('AluACA', a.gtag) in res and res[('AluACA', a.gtag)][3] > 1

    def say(cls, key, what):
        """One sentence on one row of the table, with the verdict derived from it."""
        if (cls, key) not in res:
            return f'There are no usable {what}.'
        x, y, sh, lo, hi = res[(cls, key)]
        return (f'{what}: an estimated **{_pct(sh)}** false positives '
                f'(95% CI {_pct(lo)} - {_pct(hi)}; {x:,} IP, {y:,} input calls), '
                f'**{verdict(hi)}**.')
    o.append('\nDKC1 is the H/ACA pseudouridine synthase, so snoRNA guides are the positive '
             'control, and rRNA is their canonical substrate. '
             + say('snoRNA', 'rRNA', 'snoRNA-guided rRNA chimeras') + ' '
             + say('snoRNA', a.gtag, 'snoRNA-guided genomic chimeras') + '\n')
    if have_alu:
        o.append('**The AluACA question.** '
                 + say('AluACA', a.gtag, 'AluACA-guided genomic chimeras') + '\n')
        if alu_uncal:
            o.append(f'**The method\'s assumption fails for AluACA guides.** The lower bound of '
                     f'their genomic share is {_pct(res[("AluACA", a.gtag)][3])}: the input makes '
                     'more AluACA calls per trimmed read than the IP does. If artefacts arose at the '
                     'same per-read rate in both libraries that could not happen, so for this class '
                     'they do not -- the input library is the larger source of AluACA calls. Every AluACA share in this report is therefore uncalibrated, '
                     'and the bias can differ between strata: a stratum below 100% further down is '
                     'a lead, not a measurement.\n')
        if ('AluACA', a.gtag) in res and res[('AluACA', a.gtag)][4] >= SIGNAL:
            o.append('That is not the same as saying AluACAs do not pair with DKC1 targets. '
                     'It says that, at this input depth and with the artefacts identified so '
                     'far, any real AluACA chimeras cannot be separated from the calls the '
                     'pipeline makes in the input. The rows and strata below '
                     'test whether a subset can be.\n')

    # ---- per guide -----------------------------------------------------------
    MINC = 20
    def per_guide(cls):
        gi = U_ip[U_ip.guide_class == cls]['guide_names'].value_counts()
        gc = U_ct[U_ct.guide_class == cls]['guide_names'].value_counts()
        t = pd.DataFrame({'IP': gi, 'input': gc}).fillna(0).astype(int)
        t = t[t.IP > 0]
        st = [share_ci(r.IP, r.input, N1, N2) for r in t.itertuples()]
        t['false-positive share %'] = [round(100 * x[0], 1) for x in st]
        t['upper 95% %'] = [round(100 * x[2], 1) for x in st]
        t['_hi'] = [x[2] for x in st]
        return t.sort_values('IP', ascending=False)

    guides = {c: per_guide(c) for c in ('snoRNA', 'AluACA') if c in classes}
    if have_alu:
        alu_ip, alu_ct = U_ip[U_ip.guide_class == 'AluACA'], U_ct[U_ct.guide_class == 'AluACA']
        g = guides['AluACA']
        o.append('## AluACA guides\n')
        o.append(f'**{len(alu_ip):,}** usable AluACA calls in the IP from {len(g):,} guides, '
                 f'**{len(alu_ct):,}** in the input.\n')
        crows = {}
        for c, t in guides.items():
            x = t[t.IP >= MINC]
            crows[c] = [len(x), int((x._hi < SIGNAL).sum()), int((x._hi < 0.5).sum()),
                        int((x.input == 0).sum()),
                        f'{100 * x.IP.head(5).sum() / max(t.IP.sum(), 1):.0f}%']
        o.append(md_table(pd.DataFrame.from_dict(crows, orient='index', columns=[
            f'guides with >={MINC} IP', 'upper bound < 100%', 'upper bound < 50%',
            'zero input', 'top-5 guides\' share of calls']), 'guide class'))
        sn, al = crows.get('snoRNA'), crows['AluACA']
        o.append(f"""
Calibrated against the snoRNA guides, which should look real: {sn[1]:,} of {sn[0]:,} snoRNA guides
with at least {MINC} usable IP calls have a false-positive share whose upper bound is below 100%,
against {al[1]:,} of {al[0]:,} AluACA guides. Per-guide input counts are small, so single guides
are weakly resolved either way; with one IP and one input library there is no replication
and no multiple-testing correction across {len(g):,} guides. A concentrated class -- a few guides
carrying most calls -- is a sign that particular sequences, not pairing, generate the calls.
""" if sn else '')
        o.append('Top 25 AluACA guides by usable IP calls:\n')
        o.append(md_table(g.drop(columns='_hi').head(25), 'AluACA guide'))
        o.append('')

        # ---- targets ---------------------------------------------------------
        o.append('## What the AluACA guides pair with\n')
        tc = pd.DataFrame({'IP': alu_ip['target_class'].value_counts(),
                           'input': alu_ct['target_class'].value_counts()}).fillna(0).astype(int)
        o.append(md_table(tc, 'target class (usable calls)'))
        o.append('')
        g_ip = alu_ip[alu_ip.target_class == a.gtag]
        g_ct = alu_ct[alu_ct.target_class == a.gtag]
        if not g_ip.empty:
            o.append('### Genomic arm biotypes (IP, usable)\n')
            o.append(md_table(g_ip['gene_type'].value_counts().head(12).to_frame('chimeras'), 'gene_type'))
            o.append('')
            mrna, mrna_ct = g_ip[_pc(g_ip)], g_ct[_pc(g_ct)]
            o.append('## AluACA-mRNA chimeras\n')
            o.append(f'**{len(mrna):,}** usable IP calls have a protein_coding genomic arm '
                     f'({len(mrna_ct):,} in input).\n')
            if not mrna.empty:
                feat = pd.DataFrame({'IP': mrna['feature'].value_counts(),
                                     'input': mrna_ct['feature'].value_counts()}).fillna(0).astype(int)
                o.append(md_table(feat, 'feature'))
                o.append('\nAn intronic hit is as easily explained by co-transcriptional proximity in '
                         'the host pre-mRNA as by a guide-target duplex, so intronic and exonic '
                         'counts should not be pooled.\n')

        # ---- strata ----------------------------------------------------------
        o.append('## False-positive estimate by stratum\n')
        rep_ = lambda d: d.target_in_repeat
        strata = [
            ('genomic arm, all', lambda d: d.index.notna()),
            ('arm inside a repeat', rep_),
            ('arm outside any repeat', lambda d: ~rep_(d)),
            ('outside repeat, protein_coding', lambda d: ~rep_(d) & _pc(d)),
            ('outside repeat, protein_coding, exonic', lambda d: ~rep_(d) & _pc(d) & (d.feature == 'exonic')),
        ]
        strat = {}
        for cls in ('AluACA', 'snoRNA'):
            I = U_ip[(U_ip.guide_class == cls) & (U_ip.target_class == a.gtag)]
            C = U_ct[(U_ct.guide_class == cls) & (U_ct.target_class == a.gtag)]
            rws = {}
            for lbl, fn in strata:
                x, y = int(fn(I).sum()), int(fn(C).sum())
                rws[lbl] = fp_cells(x, y, N1, N2)
                strat[(cls, lbl)] = (x, y) + share_ci(x, y, N1, N2)
            o.append(f'**{cls} guides, usable genomic calls**\n')
            o.append(md_table(pd.DataFrame.from_dict(rws, orient='index', columns=FP_COLS), 'stratum'))
            o.append('')
        sig = [lbl for lbl, _ in strata if strat[('AluACA', lbl)][4] < SIGNAL]
        if sig:
            o.append('AluACA strata whose false-positive share has an upper bound below 100%: '
                     + '; '.join(f'*{l}* ({_pct(strat[("AluACA", l)][2])}, upper '
                                 f'{_pct(strat[("AluACA", l)][4])}, {strat[("AluACA", l)][0]:,} IP calls)'
                                 for l in sig)
                     + '. These are the only places AluACA calls are IP-specific. IP-specific '
                       'includes chimeras formed after lysis, so each still needs the duplex test before '
                       'it is read as pairing.'
                     + (' Given the failed assumption above, treat them as leads: the share is '
                        'not calibrated for AluACA guides.' if alu_uncal else '') + '\n')
        else:
            o.append('**No AluACA stratum has a false-positive share whose upper bound is below '
                     '100%.** Restricting to exonic protein-coding targets outside repeats -- '
                     'the stratum a guide model predicts -- does not separate AluACA calls from '
                     'artefact either.\n')
        _EX = 'outside repeat, protein_coding, exonic'
        ex_ip = U_ip[(U_ip.guide_class == 'AluACA') & (U_ip.target_class == a.gtag) &
                     ~U_ip.target_in_repeat & _pc(U_ip) & (U_ip.feature == 'exonic')]
        own = int(ex_ip.target_in_source_locus.sum())
        o.append(f"""
**This measures trans pairing only.** {len(ex_ip) - own:,} of {len(ex_ip):,} usable exonic calls pair a guide
with an mRNA outside the guide's own locus. That is not evidence against cis action:
genome masking removes cis geometry by construction, and so does the contiguity flag -- a
guide ligated to its own host pre-mRNA yields a read that aligns contiguously and is
dropped. Testing a co-transcriptional model would need a different design.
""")

        # ---- orientation -------------------------------------------------------
        ori = alu_orientation(U_ip, U_ct, a.gtag, a.rmsk, a.bedtools,
                              os.path.join(a.resdir, 'orient_work'), N1, N2)
        if ori:
            o.append('## Alu targets, split by orientation\n')
            o.append('An Alu-derived guide can only base-pair with an *antisense* Alu; a '
                     'same-orientation Alu shares its sequence. Usable genomic calls only.\n')
            orows = {}
            for cls in ('AluACA', 'snoRNA'):
                for orient, note in (('sense', 'cannot base-pair'), ('antisense', 'can base-pair')):
                    orows[f'{cls} -> {orient} Alu ({note})'] = fp_cells(
                        ori[(cls, orient, 'ip')], ori[(cls, orient, 'ctrl')], N1, N2)
            o.append(md_table(pd.DataFrame.from_dict(orows, orient='index', columns=FP_COLS), 'stratum'))
            aa, as_ = ori[('AluACA', 'antisense', 'ip')], ori[('AluACA', 'sense', 'ip')]
            sa, ss = ori[('snoRNA', 'antisense', 'ip')], ori[('snoRNA', 'sense', 'ip')]
            if as_ and ss:
                _, pv = fisher_exact([[aa, as_], [sa, ss]])
                hi_a = share_ci(aa, ori[('AluACA', 'antisense', 'ctrl')], N1, N2)[2]
                o.append(f"""
Antisense-Alu AluACA calls, the pairable class: **{verdict(hi_a)}**. Within the IP,
antisense:sense is {aa / as_:.2f} for AluACA guides against {sa / ss:.2f} for snoRNA guides
(Fisher p = {pv:.1e}). snoRNA guides have no Alu complementarity, so theirs is the baseline
availability of antisense Alus; this comparison does not use the input.
""")

    # ---- chrM ----------------------------------------------------------------
    # Dyskerin is nuclear, so a guide:mitochondrial-RNA duplex is not physically available
    # and every usable chrM call is an artefact. The input measures the pipeline and
    # library-prep artefacts; chrM calls it does NOT explain are IP-specific artefacts --
    # chimeras formed after lysis, the kind the input cannot see.
    o.append('## chrM: a floor for chimeras formed after lysis\n')
    MITO = ['chrM', 'MT', 'chrMT']
    mrows, mres = {}, {}
    for cls in [c for c in ('snoRNA', 'AluACA') if c in classes]:
        I = U_ip[(U_ip.guide_class == cls) & (U_ip.target_class == a.gtag)]
        C = U_ct[(U_ct.guide_class == cls) & (U_ct.target_class == a.gtag)]
        mi, mc = I.reference_target.astype(str).isin(MITO), C.reference_target.astype(str).isin(MITO)
        x, y = int(mi.sum()), int(mc.sum())
        lo_, hi_ = (beta.ppf(0.025, x, len(I) - x + 1) if x else 0.0,
                    beta.ppf(0.975, x + 1, len(I) - x) if len(I) > x else 1.0)
        sh = share_ci(x, y, N1, N2)
        excess = max(0, x - y * N1 / N2)
        mrows[cls] = [len(I), x, f'{100 * x / max(len(I), 1):.2f}% ({100 * lo_:.2f}% - {100 * hi_:.2f}%)',
                      y, f'{_pct(sh[0])} ({_pct(sh[1])} - {_pct(sh[2])})' if x else '',
                      f'{excess:,.0f}', f'{100 * excess / max(len(I), 1):.2f}%']
        mres[cls] = (len(I), x, y, sh, excess)
    o.append(md_table(pd.DataFrame.from_dict(mrows, orient='index', columns=[
        'usable genomic IP', 'on chrM', 'chrM share (95% CI)', 'chrM in input',
        'false-positive share of chrM calls', 'chrM calls input cannot explain',
        'as share of usable genomic IP']), 'guide class'))
    sent = []
    for cls, (n, x, y, sh, excess) in mres.items():
        if x:
            sent.append(f'{cls}: {x:,} usable chrM calls, of which the input accounts for an estimated '
                        f'{_pct(sh[0])} (upper {_pct(sh[2])}), leaving an estimated {excess:,.0f} '
                        f'({100 * excess / max(n, 1):.2f}% of usable genomic IP calls) as IP-specific artefacts')
    o.append(f"""
A guide cannot pair with a mitochondrial RNA in vivo, so every chrM call is an artefact;
the ones the input does not explain are IP-specific artefacts, formed after lysis --
whether by crowding on the beads or by complexes meeting new RNA, which the chimeric
eCLIP method paper could not separate. {'; '.join(sent)}.

This is a floor, not an estimate of all post-lysis artefacts: it counts only partners that
happen to be mitochondrial, and mitochondrial RNAs are a fraction of what is available
for random joining. The nuclear equivalent cannot be separated from real pairing by
counting reads. Whether a guide can form the >=8 bp bipartite duplex around a target
uridine can, and that is the test these calls need next.
""")

    # ---- target list -----------------------------------------------------------
    if have_alu:
        def stratum(d):
            return d[(d.guide_class == 'AluACA') & (d.target_class == a.gtag) &
                     ~d.target_in_repeat & _pc(d) & (d.feature == 'exonic')]
        TI, TC = stratum(U_ip), stratum(U_ct)
        x_, y_, sh_, lo_s, hi_s = strat[('AluACA', _EX)]
        o.append('## AluACA-mRNA candidate list\n')
        o.append(f'Usable exonic, protein-coding calls outside any repeat: {len(TI):,} in the IP, '
                 f'{len(TC):,} in the input, estimated false-positive share {_pct(sh_)} '
                 f'(upper {_pct(hi_s)}): **{verdict(hi_s)}**.\n')

        def summarise(keys, fname, label, topn):
            t = pd.DataFrame({'IP': TI.groupby(keys).size(),
                              'input': TC.groupby(keys).size()}).fillna(0).astype(int)
            t = t[t.IP > 0]
            st = [share_ci(r.IP, r.input, N1, N2) for r in t.itertuples()]
            t['false-positive share %'] = [round(100 * v[0], 1) for v in st]
            t['upper 95% %'] = [round(100 * v[2], 1) for v in st]
            t = t.sort_values('IP', ascending=False)
            t.to_csv(os.path.join(a.resdir, fname), sep='\t')
            # Per-gene counts are too thin for a share to mean anything on screen (a gene
            # with 2 calls and none in input has an upper bound in the thousands of %),
            # so the rendered table shows counts; the TSV keeps the shares.
            show = t[['IP', 'input']].copy()
            show.index = [' / '.join(map(str, k)) if isinstance(k, tuple) else str(k) for k in show.index]
            o.append(f'### {label}\n')
            o.append(f'{len(t):,} in total; full list in `{os.path.join(reldir, fname)}`. '
                     f'Top {topn} by IP count:\n')
            o.append(md_table(show.head(topn), ' / '.join(keys)))
            o.append('')
            return t

        genes = summarise(['gene_name'], f'{a.ip}.AluACA_mRNA_targets_by_gene.tsv', 'By target gene', 30)
        pairs = summarise(['guide_names', 'gene_name'], f'{a.ip}.AluACA_mRNA_targets_by_pair.tsv',
                          'By guide-target pair', 30)
        if len(genes):
            gn = genes.index.to_series().astype(str)
            mt, rp = gn.str.contains(r'\bMT-', regex=True), gn.str.match(r'^(RPS|RPL)\d')
            n_single = int((genes.IP == 1).sum())
            o.append(f"""
**A candidate list, not identified targets.** {int(mt.sum())} mitochondrial and {int(rp.sum())} ribosomal-protein
genes carry {100 * genes[mt].IP.sum() / genes.IP.sum():.1f}% and {100 * genes[rp].IP.sum() / genes.IP.sum():.1f}% of the calls; mitochondrial ones cannot be real
(see chrM), so they mark how much of the list tracks abundance. {n_single:,} of {len(genes):,} genes
({100 * n_single / len(genes):.0f}%) rest on a single call and {('only ' + str(int((genes.IP >= 5).sum()))) if (genes.IP >= 5).any() else 'none'} have five or more, so no
per-gene share is resolved; the stratum's share above is the only aggregate statement.
Names joined by `|` are ambiguous calls, not composites.
""")
        o.append(f'{len(genes):,} distinct mRNAs and {len(pairs):,} distinct guide-target pairs.\n')

        # ---- scrutiny summary -------------------------------------------------
        o.append('## How much of the AluACA genomic set survives\n')
        G = ip[(ip.guide_class == 'AluACA') & (ip.target_class == a.gtag)]
        GU = G[G.usable]
        gtot = len(G)
        if gtot:
            o.append(f'Of {gtot:,} AluACA IP calls with a genomic arm:\n')
            o.append(f'- **{int((~G.usable).sum()):,} ({100 * (~G.usable).mean():.1f}%)** are removed '
                     'by the artefact flags (simple-repeat guide, contiguous somewhere in the '
                     'genome, too many loci, not re-aligned, or guide within 2 kb).')
            o.append(f'- of the {len(GU):,} usable, **{int(GU.target_in_repeat.sum()):,}** have the '
                     f'arm inside an annotated repeat and **{int(GU.target_in_source_locus.sum()):,}** '
                     'land on a guide locus on the same strand. Both are flagged, not removed.')
            clean = GU[~GU.target_in_repeat & ~GU.target_in_source_locus]
            ce = clean[_pc(clean) & (clean.feature == 'exonic')]
            o.append(f'\nUsable, outside repeats and off the guide loci: **{len(clean):,}**, of '
                     f'which **{len(ce):,}** are exonic protein-coding -- the most conservative '
                     'AluACA-mRNA set. Its false-positive share is in the stratum table above.\n')

    # ---- sanity check against the published run ----------------------------
    if os.path.exists(a.published):
        with open(a.published) as fh:
            pub = sum(1 for _ in fh) - 1
        ours = int((ip.target_class == a.gtag).sum())
        o.append('\n## Cross-check against the published hg19 run\n')
        o.append(f'The published output for this same sample '
                 f'(`{os.path.basename(a.published)}`, hg19, plain 1951-record snoRNA source) '
                 f'holds **{pub:,}** genomic chimeras. This rerun ({a.gtag}, merged '
                 f'2701-record source) gives **{ours:,}**.\n')
        # On an hg19 arm the build term drops out, which is the whole point of
        # running one: what is left is the catalogue and the masking index.
        same_build = a.gtag == 'hg19'
        differs = ('a source catalogue with 765 extra sequences competing for the same '
                   'reads, and a repeat-masking index built from public RepBase rather '
                   'than the lab\'s own'
                   if same_build else
                   'different genome build, a source catalogue with 765 extra sequences '
                   'competing for the same reads, and a repeat-masking index built from '
                   'public RepBase rather than the lab\'s own')
        o.append(f'\nExact agreement is not expected: {differs}. Preprocessing, however, '
                 'was validated to reproduce 4983/5000 (99.66%) of the published trimmed '
                 'sequences byte-for-byte, so any divergence is downstream of read '
                 'handling.\n')
        if same_build:
            o.append('\nBecause this arm is on the published build, the genome-build term '
                     'is removed from the comparison: the residual here is attributable to '
                     'the catalogue and the masking index alone. See the matching '
                     '`*.vs_published.txt`, which scores the overlap by read name rather '
                     'than by count.\n')

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
