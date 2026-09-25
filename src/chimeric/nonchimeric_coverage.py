#!/usr/bin/env python3
"""Per-locus-class coverage of the NON-chimeric reads, IP against input.

The chimeric pipeline keeps only reads that fail to align end-to-end to the repeat
consensus and then to the genome; everything that aligns is discarded at the masking
stage as "not a chimera". For the canonical arm that discard is 52% of the IP library and
35% of the input, and it is a standard eCLIP experiment: the reads say which RNAs the bait
was crosslinked to, with no ligation step involved. This script reads the masking BAM the
pipeline leaves behind and answers the question the chimera tables cannot:

    are reads over the AluACA loci enriched in the IP over the input,
    the way reads over the snoRNA loci are?

That is a binding question, not a pairing question. It does not depend on the chimeric
ligation, on the guide/target arm split, or on --outFilterMultimapNmax 1.

It also emits the library-size table, which the chimeric result needs as a normaliser:
the two libraries neither trim nor map alike (79% vs 44% surviving trimming, 50.6% vs
33.3% mapping to the genome), so "per trimmed read" is not a like-for-like denominator
and the false-positive shares in RESULTS.md move materially with the choice.

Three things are reported side by side rather than chosen, because the answer depends on
each and a silent choice would hide that:

  multimapper policy   unique (NH==1) | primary | fractional (1/NH)
                       The masking step ran --outFilterMultimapNmax 100, so an Alu read
                       typically has many alignments. `unique` undercounts every repeat
                       family; `fractional` is what the published repeatquant does in
                       spirit; `primary` is arbitrary among ties and is shown only so a
                       reader can see it is arbitrary.
  orientation          same-strand | opposite-strand relative to the locus.
                       Not assumed: the snoRNA loci are the internal calibrator, since
                       they must be overwhelmingly one orientation. The report names
                       which one the data chose.
  duplicates           raw | UMI-deduplicated by (UMI, chrom, pos, strand).
                       The pipeline extracts UMIs into the read name but never dedups, so
                       raw is what the chimera counts in RESULTS.md are comparable to.

Usage:
    python3 src/chimeric/nonchimeric_coverage.py \\
        --resdir results/chimeric_chr25/arm3_hg19_merged \\
        --ip SRR30692552 --input SRR30692553 --gtag hg19 \\
        --loci-bed ref/chimeric/guide_loci.hg19.bed \\
        --out results/chimeric_chr25/nonchimeric_coverage.tsv

The BAM it needs is $resdir/<uid>/<uid>.genome.map/Aligned.out.bam, which sno-chimeras.py
keeps because run_chimeras.sh passes --keep. If a run was fetched from another machine
without it, re-run just the masking step rather than the whole pipeline.
"""

import argparse
import bisect
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import find_input, out_path, require          # noqa: E402

CONF = 0.95
POLICIES = ('unique', 'primary', 'fractional')
# Walk-back bound for the interval lookup. The longest union interval is 6,859 nt; this
# leaves two orders of magnitude of headroom so the scan cannot miss a long locus.
MAX_LOCUS = 1_000_000


# --- statistics ---------------------------------------------------------------

def rate_ratio_ci(x, y, n1, n2, conf=CONF):
    """IP/input rate ratio (x/n1)/(y/n2) with an exact conditional-binomial interval.

    x ~ Pois(l1*n1), y ~ Pois(l2*n2); conditional on x+y, x ~ Binomial(n, p) with
    ratio = l1/l2 = (p/(1-p)) * (n2/n1), so a Clopper-Pearson interval on p maps
    monotonically onto the ratio. Same construction as make_report.share_ci, inverted:
    there the statistic is the input's share of the IP, here it is the enrichment.

    Counts may be fractional under the 1/NH policy, which the Poisson model does not
    strictly allow; the interval is then approximate and is labelled as such.
    """
    from scipy.stats import beta
    if y == 0 and x == 0:
        return float('nan'), float('nan'), float('nan')
    n = x + y
    a = (1 - conf) / 2
    p_lo = 0.0 if x == 0 else beta.ppf(a, x, n - x + 1)
    p_hi = 1.0 if y == 0 else beta.ppf(1 - a, x + 1, n - x)
    f = lambda q: float('inf') if q >= 1 else (q / (1 - q)) * (n2 / n1)
    point = float('inf') if y == 0 else (x / n1) / (y / n2)
    return point, f(p_lo), f(p_hi)


# --- inputs -------------------------------------------------------------------

def read_counts(outdir, uid):
    """Per-stage read counts from the logs, and the non-chimeric totals derived from them.

    Mirrors make_report.read_counts and extends it: that function reports what SURVIVES
    masking (the chimera candidates), this one also reports what masking REMOVED, which is
    the non-chimeric fraction this script is about.
    """
    n = {}
    cut = os.path.join(outdir, f'{uid}.cut.adapt.log')
    if os.path.exists(cut):
        with open(cut) as fh:
            rows = [l.split('\t') for l in fh if l.startswith('OK')]
        if rows:
            n['raw'] = int(rows[0][1])
            n['trimmed'] = int(rows[-1][6])
    for stage, fn in (('repeat', f'{uid}.mask.repeat.map.log'),
                      ('genome', f'{uid}.mask.genome.map.log')):
        p = os.path.join(outdir, fn)
        if not os.path.exists(p):
            continue
        txt = open(p).read()

        def grab(pat):
            m = re.search(pat + r'\s*\|\s*([0-9.]+)', txt)
            return float(m.group(1)) if m else None

        total = grab(r'Number of input reads')
        uniq = grab(r'Uniquely mapped reads number') or 0
        multi = grab(r'Number of reads mapped to multiple loci') or 0
        if total is None:
            continue
        n[f'{stage}_in'] = int(total)
        n[f'{stage}_mapped'] = int(uniq + multi)
        n[f'{stage}_unique'] = int(uniq)
    if 'genome_in' in n:
        n['candidates'] = n['genome_in'] - n['genome_mapped']
        n['nonchimeric'] = n.get('repeat_mapped', 0) + n['genome_mapped']
    return n


def _load_union(path):
    """name -> source ('naprnadb_only', 'jady_aluaca', ...), or name -> None.

    Read by COLUMN NAME, never by position. The two union files disagree and step 09
    changes both: AluACA_union_nr.bed is BED6 before step 09 and BED6+repeat_family after
    it, so its 7th column is the repeat family, not the source -- reading position 7 would
    silently classify every NapRNAdb-only locus as Jady-backed. Only the .tsv carries a
    header, so a header is required to read `source`; without one the names are still used
    for AluACA membership but the split is refused rather than guessed.
    """
    if not path or not os.path.exists(path):
        return {}, False
    with open(path) as fh:
        first = fh.readline().rstrip('\n').split('\t')
        header = {c: i for i, c in enumerate(first)} if first and first[0] == 'chrom' else None
        if header is None:
            fh.seek(0)
        i_name = header['name'] if header and 'name' in header else 3
        i_src = header.get('source') if header else None
        out = {}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) <= i_name:
                continue
            out[f[i_name]] = f[i_src] if i_src is not None and len(f) > i_src else None
    return out, i_src is not None


def load_loci(loci_bed, union_bed, snodb_tsv):
    """Locus intervals plus a class per locus.

    AluACA membership is decided by exact name match against the union BED, never by the
    `.id3xxx` pattern -- annotate_chimeras.load_alu_names documents why. The union's own
    `source` column splits Jady-backed loci from NapRNAdb-only ones, which matters because
    those two halves differ in every other respect. snoDB box_type refines the rest into
    H/ACA and C/D when the step-10 TSV is available, since only H/ACA is a DKC1 substrate
    class and C/D is therefore a second internal negative control.
    """
    alu_src, have_source = _load_union(union_bed)

    box = {}
    if snodb_tsv and os.path.exists(snodb_tsv):
        with open(snodb_tsv) as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) >= 8 and f[0] != 'chrom':
                    box[f[3]] = f[7]

    def classify(name):
        if name in alu_src:
            src = alu_src[name]
            if src is None:
                return 'AluACA_unsplit'
            return 'AluACA_naprnadb_only' if src == 'naprnadb_only' else 'AluACA_jady'
        b = box.get(name, '')
        if b in ('H/ACA', 'AluACA', 'Alu-ACA'):
            return 'snoRNA_HACA'
        if b in ('C/D', 'SNORD-like'):
            return 'snoRNA_CD'
        if b == 'scaRNA':
            return 'scaRNA'
        return 'other'

    by_chrom = defaultdict(list)
    n = 0
    with open(loci_bed) as fh:
        for line in fh:
            if line.startswith(('#', 'track', 'chrom\t')):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 6:
                continue
            by_chrom[f[0]].append((int(f[1]), int(f[2]), f[3], f[5], classify(f[3])))
            n += 1
    if not n:
        sys.exit(f'no usable BED6 records in {loci_bed}')
    index = {}
    for c, rows in by_chrom.items():
        rows.sort()
        index[c] = (rows, [r[0] for r in rows])
    return index, n


# --- counting -----------------------------------------------------------------

def count_bam(bam, index, dedup=True, min_overlap=10):
    """One pass over the masking BAM, accumulating per locus.

    Single pass with a bisect lookup rather than a fetch per region: it needs no sorted,
    indexed BAM (STAR writes the masking BAM unsorted) and no temp files. A read
    overlapping two loci is counted in both; the per-locus TSV makes that auditable.

    Deduplication is applied only to reads that touch a locus, so the set stays small.
    The library totals come from the masking logs and are NOT deduplicated -- which is the
    right comparison for RESULTS.md, whose chimera counts are not deduplicated either.
    """
    import pysam
    per_locus = defaultdict(lambda: defaultdict(float))
    seen = set()
    n_aln = n_hit = 0
    with pysam.AlignmentFile(bam, 'rb', check_sq=False) as fh:
        for a in fh:
            if a.is_unmapped:
                continue
            n_aln += 1
            chrom = a.reference_name
            if chrom not in index:
                continue
            rows, starts = index[chrom]
            # bisect_left, not _right: a locus starting exactly at reference_end does not
            # overlap a half-open interval.
            i = bisect.bisect_left(starts, a.reference_end)
            aln_strand = '-' if a.is_reverse else '+'
            nh = a.get_tag('NH') if a.has_tag('NH') else 1
            umi = a.query_name.rsplit('_', 1)[-1]
            for start, end, name, strand, cls in reversed(rows[:i]):
                # rows are start-sorted and walked backwards, so once a start is far
                # enough upstream that no locus of any plausible length could reach the
                # read, nothing earlier can either.
                if start + MAX_LOCUS < a.reference_start:
                    break
                ov = min(end, a.reference_end) - max(start, a.reference_start)
                if ov < min(min_overlap, end - start):
                    continue
                n_hit += 1
                orient = 'same' if aln_strand == strand else 'opposite'
                if dedup:
                    key = (umi, chrom, a.reference_start, aln_strand, name)
                    if key in seen:
                        continue
                    seen.add(key)
                d = per_locus[(name, cls, orient)]
                d['fractional'] += 1.0 / nh
                if nh == 1:
                    d['unique'] += 1
                if not a.is_secondary and not a.is_supplementary:
                    d['primary'] += 1
    return per_locus, n_aln, n_hit


# --- reporting ----------------------------------------------------------------

def fmt(v):
    if v != v:
        return 'n/a'
    return 'inf' if v == float('inf') else f'{v:.2f}'


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--resdir', default=os.environ.get('OUT', '') or None,
                   help='arm directory holding <uid>/ per sample')
    p.add_argument('--ip', default='SRR30692552')
    p.add_argument('--input', dest='inp', default='SRR30692553')
    p.add_argument('--ip-bam', default=None, help='override the IP masking BAM path')
    p.add_argument('--input-bam', default=None, help='override the input masking BAM path')
    p.add_argument('--gtag', default=os.environ.get('SPECIES', 'hg19'),
                   help='genome tag; picks the default loci BED')
    p.add_argument('--loci-bed', default=None,
                   help='BED6 of guide loci [ref/chimeric/guide_loci.<gtag>.bed]')
    p.add_argument('--union-bed', default=None,
                   help='AluACA union table, for AluACA membership and the jady/NapRNAdb '
                        'split. Use the .tsv: only it names its columns, and the .bed\'s '
                        '7th column is repeat_family after step 09, not source. '
                        '[AluACA_union_nr.tsv]')
    p.add_argument('--snodb', default=None,
                   help='snoDB_with_AluACA_union.tsv, for box_type classes (optional)')
    p.add_argument('--no-dedup', action='store_true',
                   help='skip UMI deduplication of locus-overlapping reads')
    p.add_argument('--min-overlap', type=int, default=10,
                   help='minimum bp a read must share with a locus, capped at the locus '
                        'length so short loci are not excluded [10]')
    p.add_argument('--out', default=None, help='per-locus TSV [<resdir>/nonchimeric_coverage.tsv]')
    a = p.parse_args()

    if not a.resdir:
        sys.exit('--resdir is required (the arm directory holding <uid>/)')
    loci_bed = a.loci_bed or find_input(f'guide_loci.{a.gtag}.bed')
    union_bed = (a.union_bed or os.environ.get('ALU_TSV')
                 or find_input('AluACA_union_nr.tsv'))
    snodb = a.snodb or find_input('snoDB_with_AluACA_union.tsv')
    bams = {}
    for lab, uid, override in (('IP', a.ip, a.ip_bam), ('input', a.inp, a.input_bam)):
        bams[lab] = override or os.path.join(a.resdir, uid, f'{uid}.genome.map',
                                            'Aligned.out.bam')
    require([('loci BED (--loci-bed)', loci_bed)]
            + [(f'{lab} masking BAM (--{lab.lower()}-bam); sno-chimeras.py keeps it '
                f'under <uid>.genome.map/ when run with --keep', b)
               for lab, b in bams.items()])

    index, n_loci = load_loci(loci_bed, union_bed, snodb)
    print(f'loci: {n_loci} from {loci_bed}')
    if not os.path.exists(union_bed):
        print(f'  !! no union table at {union_bed}: AluACA loci cannot be identified, '
              f'every locus falls in "other". Pass --union-bed.')
    elif not _load_union(union_bed)[1]:
        print(f'  !! {union_bed} has no `source` column, so Jady-backed and NapRNAdb-only '
              f'loci\n     cannot be separated; they are reported together as '
              f'AluACA_unsplit. Pass the .tsv.')

    lib, counts = {}, {}
    for lab, uid in (('IP', a.ip), ('input', a.inp)):
        lib[lab] = read_counts(os.path.join(a.resdir, uid), uid)
        if 'nonchimeric' not in lib[lab]:
            sys.exit(f'{lab}: could not read masking logs under '
                     f'{os.path.join(a.resdir, uid)}; they are what gives the library size')
        per_locus, n_aln, n_hit = count_bam(bams[lab], index, dedup=not a.no_dedup,
                                            min_overlap=a.min_overlap)
        counts[lab] = per_locus
        print(f'{lab}: {n_aln:,} alignments in the masking BAM, {n_hit:,} locus overlaps')

    # ---- library-size table: the normaliser the chimeric result needs
    print('\n=== library sizes (from the masking logs) ===')
    keys = [('raw', 'raw reads'), ('trimmed', 'after trimming'),
            ('repeat_mapped', 'repeat-mapped (masked)'),
            ('genome_mapped', 'genome-mapped end-to-end'),
            ('genome_unique', '  of which uniquely'),
            ('nonchimeric', 'NON-CHIMERIC total'),
            ('candidates', 'chimera candidates')]
    print(f'  {"":26s} {"IP":>12s} {"input":>12s}   {"IP%":>6s} {"input%":>7s}')
    for k, lab in keys:
        vi, vn = lib['IP'].get(k), lib['input'].get(k)
        if vi is None or vn is None:
            continue
        ti, tn = lib['IP'].get('trimmed', 1), lib['input'].get('trimmed', 1)
        print(f'  {lab:26s} {vi:>12,} {vn:>12,}   {100*vi/ti:5.1f}% {100*vn/tn:6.1f}%')
    print('\n  Candidate denominators for the IP/input comparison in RESULTS.md:')
    for k, lab in (('trimmed', 'trimmed reads (what make_report uses)'),
                   ('candidates', 'post-masking chimera candidates'),
                   ('genome_mapped', 'non-chimeric genome-mapped reads')):
        if k in lib['IP']:
            print(f'    {lab:40s} IP {lib["IP"][k]:>10,}  input {lib["input"][k]:>10,}')

    # ---- orientation, calibrated on the snoRNA loci rather than assumed
    print('\n=== orientation convention, calibrated on the snoRNA loci ===')
    orient_tot = {}
    for o in ('same', 'opposite'):
        orient_tot[o] = sum(v['fractional'] for (nm, c, oo), v in counts['IP'].items()
                            if oo == o and c.startswith('snoRNA'))
    tot = sum(orient_tot.values())
    if tot:
        best = max(orient_tot, key=orient_tot.get)
        for o in ('same', 'opposite'):
            print(f'  snoRNA loci, {o:9s}-strand reads: {orient_tot[o]:12,.0f}  '
                  f'({100*orient_tot[o]/tot:5.1f}%)')
        print(f'  -> using {best!r}-strand as the RNA orientation for this library')
    else:
        best = 'same'
        print('  no snoRNA-class reads; falling back to same-strand (UNCALIBRATED)')

    # ---- per-class enrichment
    n1k, n2k = 'genome_mapped', 'genome_mapped'
    n1, n2 = lib['IP'][n1k], lib['input'][n2k]
    print(f'\n=== per-class IP/input enrichment ({best}-strand, per million '
          f'non-chimeric genome-mapped reads) ===')
    classes = sorted({c for d in counts.values() for (_, c, _) in d})
    for pol in POLICIES:
        print(f'\n  -- multimapper policy: {pol}')
        print(f'     {"class":24s} {"IP":>11s} {"input":>10s} {"IP/M":>9s} {"input/M":>9s}'
              f' {"enrichment (95% CI)":>26s}')
        for cls in classes:
            x = sum(v[pol] for (nm, c, o), v in counts['IP'].items()
                    if c == cls and o == best)
            y = sum(v[pol] for (nm, c, o), v in counts['input'].items()
                    if c == cls and o == best)
            if x == 0 and y == 0:
                continue
            r, lo, hi = rate_ratio_ci(x, y, n1, n2)
            sym = 'inf' if r == float('inf') else f'{fmt(r)}x'
            ci = f'{sym} ({fmt(lo)} - {fmt(hi)})'
            print(f'     {cls:24s} {x:>11,.0f} {y:>10,.0f} {1e6*x/n1:>9.1f} '
                  f'{1e6*y/n2:>9.1f} {ci:>26s}')
        if pol == 'fractional':
            print('     (fractional counts are not integer Poisson draws, so these '
                  'intervals are approximate)')

    print('\nRead this against the snoRNA_HACA row, not in isolation: that is the class '
          '\n DKC1 must be enriched on. If it is not enriched here, nothing else in this '
          '\n table is interpretable. snoRNA_CD is a second negative control -- C/D box '
          '\n snoRNAs are not DKC1 substrates.')

    # ---- per-locus TSV
    out = a.out or os.path.join(a.resdir, 'nonchimeric_coverage.tsv')
    names = {}
    for lab in counts:
        for (nm, cls, o), v in counts[lab].items():
            names.setdefault((nm, cls), {})[(lab, o)] = v
    with open(out, 'w') as fh:
        cols = ['name', 'class'] + [f'{lab}_{o}_{pol}' for lab in ('IP', 'input')
                                    for o in ('same', 'opposite') for pol in POLICIES]
        fh.write('\t'.join(cols) + '\n')
        for (nm, cls), d in sorted(names.items()):
            row = [nm, cls]
            for lab in ('IP', 'input'):
                for o in ('same', 'opposite'):
                    v = d.get((lab, o), {})
                    row += [f'{v.get(pol, 0):.3f}' for pol in POLICIES]
            fh.write('\t'.join(row) + '\n')
    print(f'\nper-locus counts -> {out}  ({len(names)} loci with any read)')


if __name__ == '__main__':
    main()
