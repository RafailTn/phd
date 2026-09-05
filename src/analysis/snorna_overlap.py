#!/usr/bin/env python3
"""
Sequence overlap between two FASTA sets - built for the AluACA union against a
snoRNA catalogue (snoRNA.txt.fa), but nothing here is specific to either.

Four passes, cheapest first:

  exact        whole-string equality.  A dict lookup, so a single nucleotide of
               difference in the interval boundaries is a miss.  Understates
               agreement badly: two catalogues calling the same locus with
               113 vs 114 nt intervals are not "identical" by this test.
  containment  one sequence is a substring of the other.  Still exact character
               matching, but tolerant of different interval extents.  This is
               the pass that finds same-locus entries.
  revcomp      exact equality after reverse-complementing.  Catches a set
               deposited on the opposite strand.
  align        Smith-Waterman (EMBOSS `water`), each surviving B entry against
               every A.  The only pass that sees substitutions.  Optional -
               needs `water` on PATH, or --water /path/to/water.

The align pass is prefiltered: a B entry is aligned only if it shares an exact
k-mer (default k=20) with some A sequence.  Aligning all of B is one `water`
invocation per record and takes hours.  Two sequences above ~90% identity over
100 nt are near-certain to share an exact 20-mer, so the filter is safe at the
identities worth reporting; it does get lossy in the 75-85% band, where the
hits are Alu-ancestry background rather than shared loci.  Lower --kmer to
widen it, at a cost in runtime.

Report the exact and containment counts separately.  They answer different
questions and the gap between them is boundary convention, not biology.

DO NOT score similarity with difflib.SequenceMatcher.  Its autojunk heuristic
treats any element occurring in >1% of a sequence of 200+ elements as junk; on
a 4-letter alphabet that is every base, so sequences >=200 nt silently score
0.0 and drop out of the ranking.  An earlier version of this analysis used it
and invented a bimodal identity distribution that does not exist.

Usage:
  python3 snorna_overlap.py --a AluACA_union_nr.fasta --b snoRNA.txt.fa
  python3 snorna_overlap.py --a ... --b ... --align --out overlap_report.tsv

`water` is not in deps/pixi.toml.  Either add it:
    pixi add emboss
or run this script under a throwaway env:
    pixi exec -c conda-forge -c bioconda --spec emboss -- python3 ...
"""
import argparse
import collections
import os
import re
import shutil
import subprocess
import sys
import tempfile

# paths.py lives one level up, shared with the chimeric pipeline; the repo is a
# collection of scripts rather than an installed package, so put src/ on the path.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import find_input, find_tool, out_path, require

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    """{name: sequence}. Name is the header up to the first whitespace.

    Uppercases, so a soft-masked reference and an unmasked one compare equal.
    """
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def check_alphabet(label, seqs):
    """A U/T mismatch between the two files would make every test fail
    silently, so say what is actually in them."""
    chars = collections.Counter(c for v in seqs.values() for c in v)
    odd = {c: n for c, n in chars.items() if c not in "ACGT"}
    print(f"  {label:24s} {len(seqs):5d} records, "
          f"{len(set(seqs.values())):5d} distinct sequences"
          + (f", non-ACGT: {odd}" if odd else ""))
    if "U" in chars:
        print(f"    WARNING: {label} is RNA alphabet (U); the other set is "
              f"probably DNA (T). Convert before comparing.", file=sys.stderr)


def within_file_duplicates(label, seqs):
    """Sequences held by more than one record. These collapse in a
    sequence-keyed index, so report them rather than lose them."""
    m = collections.defaultdict(list)
    for k, v in seqs.items():
        m[v].append(k)
    dup = [v for v in m.values() if len(v) > 1]
    if dup:
        print(f"  {label}: {len(dup)} sequence(s) held by >1 record")
        for names in dup:
            print(f"      {', '.join(names)}")
    return dup


def kmer_prefilter(a_seqs, b_seqs, k):
    """B entries sharing at least one exact k-mer with some A sequence."""
    index = set()
    for v in a_seqs.values():
        for i in range(len(v) - k + 1):
            index.add(v[i:i + k])
    keep = {}
    for bk, bv in b_seqs.items():
        if any(bv[i:i + k] in index for i in range(len(bv) - k + 1)):
            keep[bk] = bv
    return keep


def water_best(a_seqs, b_seqs, water, min_cov):
    """Best Smith-Waterman hit in A for each B entry.

    One `water` invocation per B entry against the whole of A. Hits covering
    less than min_cov of the shorter sequence are discarded: a 30 nt alignment
    at 100% identity between two 150 nt Alu-derived sequences means only that
    both are Alus.

    Best is highest percent identity, then longest alignment. Ranking by raw
    identical-base count instead would let a long, mediocre alignment outrank a
    short perfect one and put a contained sequence below 100% - which reads as
    a disagreement with the containment pass when there is none. The min_cov
    floor is what keeps a short spurious alignment from winning.
    """
    out = {}
    with tempfile.TemporaryDirectory() as td:
        target = os.path.join(td, "a.fa")
        with open(target, "w") as fh:
            for k, v in a_seqs.items():
                fh.write(f">{k}\n{v}\n")
        for i, (bk, bv) in enumerate(sorted(b_seqs.items()), 1):
            query = os.path.join(td, "q.fa")
            with open(query, "w") as fh:
                fh.write(f">{bk}\n{bv}\n")
            res = os.path.join(td, "w.txt")
            subprocess.run(
                [water, "-asequence", query, "-bsequence", target,
                 "-gapopen", "10", "-gapextend", "0.5",
                 "-outfile", res, "-aformat3", "pair"],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            best, cur = None, None
            for line in open(res):
                m = re.match(r"^# 2: (\S+)", line)
                if m:
                    cur = m.group(1)
                m = re.match(r"^# Identity:\s+(\d+)/(\d+)", line)
                if m:
                    ident, alen = int(m.group(1)), int(m.group(2))
                    cov = alen / min(len(bv), len(a_seqs[cur]))
                    pid = 100.0 * ident / alen
                    if cov >= min_cov and (best is None or (pid, alen) > (best[3], best[2])):
                        best = (cur, ident, alen, pid, cov)
            if best:
                out[bk] = best
            print(f"\r  aligning {i}/{len(b_seqs)} ...", end="", file=sys.stderr)
        print("", file=sys.stderr)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--a", default=find_input("AluACA_union_nr.fasta"),
                   help="FASTA A [AluACA_union_nr.fasta]")
    p.add_argument("--b", default=find_input("snoRNA.txt.fa"),
                   help="FASTA B [snoRNA.txt.fa]")
    p.add_argument("--out", default=out_path("AluACA_snoRNA_overlap.tsv"),
                   help="write the per-pair report here (TSV)")
    p.add_argument("--align", action="store_true",
                   help="also run Smith-Waterman (needs EMBOSS water)")
    p.add_argument("--water", default=None, help="path to EMBOSS water")
    p.add_argument("--kmer", type=int, default=20,
                   help="prefilter for --align: only align B entries sharing an "
                        "exact k-mer with A [20]")
    p.add_argument("--min-cov", type=float, default=0.8,
                   help="minimum alignment coverage of the shorter sequence "
                        "for an SW hit to count [0.8]")
    args = p.parse_args()

    require([("FASTA A (--a)", args.a), ("FASTA B (--b)", args.b)])

    A, B = read_fasta(args.a), read_fasta(args.b)
    print("inputs")
    check_alphabet(os.path.basename(args.a), A)
    check_alphabet(os.path.basename(args.b), B)
    within_file_duplicates(os.path.basename(args.a), A)
    within_file_duplicates(os.path.basename(args.b), B)

    # index B by sequence; a list, because two records may share a sequence
    b_index = collections.defaultdict(list)
    for k, v in B.items():
        b_index[v].append(k)

    rows = []
    exact = set()
    for ak, av in A.items():
        for bk in b_index.get(av, ()):
            rows.append((ak, bk, "exact", len(av), len(B[bk]), len(av), "100.0", "1.00"))
            exact.add((ak, bk))

    revcomp = set()
    for ak, av in A.items():
        for bk in b_index.get(rc(av), ()):
            rows.append((ak, bk, "revcomp", len(av), len(B[bk]), len(av), "100.0", "1.00"))
            revcomp.add((ak, bk))

    # containment is the O(n*m) pass; short sequences are skipped because a
    # 20 nt substring hit between two Alus carries no information
    contained = set()
    for bk, bv in B.items():
        if len(bv) < 40:
            continue
        for ak, av in A.items():
            if (ak, bk) in exact:
                continue
            if bv in av:
                kind, ov = "b_in_a", len(bv)
            elif av in bv:
                kind, ov = "a_in_b", len(av)
            else:
                continue
            rows.append((ak, bk, kind, len(av), len(bv), ov, "100.0",
                         f"{ov / min(len(av), len(bv)):.2f}"))
            contained.add((ak, bk))

    aligned = {}
    if args.align:
        water = find_tool("water", args.water)
        if not water:
            print("\n--align requested but EMBOSS `water` was not found. "
                  "Install it (pixi add emboss) or pass --water PATH.",
                  file=sys.stderr)
            sys.exit(1)
        candidates = kmer_prefilter(A, B, args.kmer)
        print(f"\nSmith-Waterman: {len(candidates)}/{len(B)} B entries share an "
              f"exact {args.kmer}-mer with A and will be aligned")
        aligned = water_best(A, candidates, water, args.min_cov)
        for bk, (ak, ident, alen, pid, cov) in aligned.items():
            if (ak, bk) in exact or (ak, bk) in contained:
                continue
            rows.append((ak, bk, "align", len(A[ak]), len(B[bk]), alen,
                         f"{pid:.1f}", f"{cov:.2f}"))

    print("\nresults")
    print(f"  identical sequences:            {len(exact)}")
    print(f"  reverse-complement identical:   {len(revcomp)}")
    print(f"  one wholly inside the other:    {len(contained)}")
    print(f"  same locus (exact+contained):   {len(exact | contained)}")
    if aligned:
        hist = collections.Counter()
        for ak, ident, alen, pid, cov in aligned.values():
            hist[min(int(pid // 5) * 5, 100)] += 1
        print(f"\n  SW best hit per B entry (coverage >= {args.min_cov:.0%} "
              f"of the shorter sequence), by identity:")
        for band in sorted(hist, reverse=True):
            print(f"    {band:3d}-{band + 4}%  {hist[band]:4d}  {'#' * min(hist[band], 60)}")

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("a_name\tb_name\trelation\ta_len\tb_len\toverlap_len\t"
                     "identity_pct\tcoverage\n")
            for r in sorted(rows, key=lambda r: (-float(r[6]), -int(r[5]))):
                fh.write("\t".join(str(x) for x in r) + "\n")
        print(f"\nwrote {len(rows)} pair(s) to {args.out}")


if __name__ == "__main__":
    main()
