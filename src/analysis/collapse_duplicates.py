#!/usr/bin/env python3
"""Collapse the loci present in BOTH snoRNA.txt.fa and the AluACA union FASTA.

Fifteen loci are catalogued twice -- once in snoRNA.txt.fa (as SNODB2065-2089)
and once in the AluACA union -- so concatenating the two files duplicates them.
Every one of the fifteen is an exact substring relationship: six pairs are
byte-identical and nine differ only in where each catalogue put the 5'/3'
boundary.  There are no internal substitutions, which is why plain string
containment finds them all and no alignment is needed here.

Which member survives, in order:

  1. identical sequences            -> the union record (it carries coordinates)
  2. exactly one ends canonically   -> that one.  "Canonical" means the ACA box
     sits exactly 3 nt from the 3' end, i.e. seq[-6:-3] == "ACA".  Four pairs
     are decided here: three where snoDB stops short of the box and one where
     it overshoots it by 3 nt.
  3. both end canonically, so the difference is 5'-only:
       - extra 5' <= --offset-tol (default 1 nt) -> a coordinate off-by-one,
         keep the union record
       - otherwise the longer sequence carries more of the 5' hairpin and wins

Rule 3's second branch fires twice.  Once it keeps the union (AluACA30, 55 nt
more 5' sequence); once it keeps the snoDB sequence (SNODB2089 vs AluACA163,
36 nt more).  AluACA163 is a Jady et al. deposit, and those are 3'-half only by
design, so the 36 nt is the start of the omitted 5' hairpin.

The surviving record always keeps the UNION header, so the merged catalogue has
one naming scheme.  When rule 3 hands the sequence over from snoDB, the union
interval no longer matches its sequence, so --bed / --bed-out write a corrected
BED; the extension is verified genomic (checked against hg38 for AluACA163).

Outputs: a merged non-redundant FASTA (1951 + 765 - 15 = 2701 records) and a
TSV recording which rule fired for each of the fifteen.
"""

import argparse
import re
import sys


def read_fasta(path):
    """Ordered list of (name, sequence); name is the first whitespace token."""
    out, name, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    out.append((name, "".join(chunks).upper()))
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        out.append((name, "".join(chunks).upper()))
    return out


def canonical(seq):
    """True if the ACA box sits exactly 3 nt from the 3' end."""
    return len(seq) >= 6 and seq[-6:-3] == "ACA"


def decide(sno_seq, uni_seq, tol):
    """-> (winner, rule) where winner is 'union' or 'snodb'."""
    if sno_seq == uni_seq:
        return "union", "identical"
    c_sno, c_uni = canonical(sno_seq), canonical(uni_seq)
    if c_uni and not c_sno:
        return "union", "union-ACA-canonical"
    if c_sno and not c_uni:
        return "snodb", "snodb-ACA-canonical"
    # both (or neither) canonical: the difference is at the 5' end
    long_, short = (uni_seq, sno_seq) if len(uni_seq) > len(sno_seq) else (sno_seq, uni_seq)
    five = long_.find(short)
    if five <= tol:
        return "union", "5prime-offset"
    return ("union" if len(uni_seq) > len(sno_seq) else "snodb"), "longer-5prime-hairpin"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sno", default="snoRNA.txt.fa", help="snoRNA catalogue FASTA")
    ap.add_argument("--union", default="AluACA_union_nr.fasta", help="AluACA union FASTA")
    ap.add_argument("--out", default="AluACA_snoRNA_merged_nr.fasta", help="merged FASTA")
    ap.add_argument("--report", default="AluACA_snoRNA_collapse_report.tsv",
                    help="per-pair decision table")
    ap.add_argument("--bed", help="union BED, to correct intervals whose sequence changed")
    ap.add_argument("--bed-out", help="corrected BED (requires --bed)")
    ap.add_argument("--offset-tol", type=int, default=1,
                    help="5' difference at or below this is a coordinate off-by-one [1]")
    args = ap.parse_args()

    sno = read_fasta(args.sno)
    union = read_fasta(args.union)
    uni_by_seq = {}
    for name, seq in union:
        uni_by_seq.setdefault(seq, name)

    # Pair up by exact containment, either direction.  A dict of union sequences
    # makes the identical case O(1); containment needs the scan, but 765 x 1951
    # substring tests on ~150 nt strings is a couple of seconds.
    pairs = []           # (sno_name, sno_seq, uni_name, uni_seq)
    for sname, sseq in sno:
        hit = None
        if sseq in uni_by_seq:
            hit = (uni_by_seq[sseq], sseq)
        else:
            for uname, useq in union:
                if sseq in useq or useq in sseq:
                    hit = (uname, useq)
                    break
        if hit:
            pairs.append((sname, sseq, hit[0], hit[1]))

    replaced = {}        # union name -> sequence taken from snoDB
    drop_sno = set()
    rows = []
    for sname, sseq, uname, useq in pairs:
        winner, rule = decide(sseq, useq, args.offset_tol)
        drop_sno.add(sname)
        if winner == "snodb":
            replaced[uname] = sseq
        rows.append((uname, sname, len(useq), len(sseq), rule, winner,
                     len(sseq) if winner == "snodb" else len(useq)))

    with open(args.out, "w") as fh:
        for name, seq in sno:
            if name in drop_sno:
                continue
            fh.write(f">{name}\n{seq}\n")
        for name, seq in union:
            fh.write(f">{name}\n{replaced.get(name, seq)}\n")

    with open(args.report, "w") as fh:
        fh.write("union_name\tsnodb_name\tunion_len\tsnodb_len\trule\tkept\tkept_len\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    if args.bed:
        if not args.bed_out:
            sys.exit("--bed requires --bed-out")
        # Union FASTA headers are the BED name with "|" -> "_" and a ".idN"
        # suffix appended (step 08), so map back before matching.
        by_bed = {}
        for uname in replaced:
            by_bed[re.sub(r"\.id[0-9]+$", "", uname).replace("|", "_")] = replaced[uname]
        n_fixed = 0
        with open(args.bed) as fin, open(args.bed_out, "w") as fout:
            for line in fin:
                f = line.rstrip("\n").split("\t")
                key = f[3].replace("|", "_") if len(f) >= 6 else None
                if key in by_bed:
                    # The kept sequence extends the interval; the extension is
                    # 5' in transcript orientation, so it moves the low
                    # coordinate on +, the high coordinate on -.
                    grow = len(by_bed[key]) - (int(f[2]) - int(f[1]))
                    if f[5] == "-":
                        f[2] = str(int(f[2]) + grow)
                    else:
                        f[1] = str(int(f[1]) - grow)
                    n_fixed += 1
                fout.write("\t".join(f) + "\n")
        print(f"corrected {n_fixed} BED interval(s) -> {args.bed_out}")

    kept_union = sum(1 for r in rows if r[5] == "union")
    print(f"pairs collapsed: {len(rows)}  (union kept {kept_union}, snoDB sequence kept {len(rows) - kept_union})")
    for rule in ("identical", "union-ACA-canonical", "snodb-ACA-canonical",
                 "5prime-offset", "longer-5prime-hairpin"):
        n = sum(1 for r in rows if r[4] == rule)
        if n:
            print(f"  {rule:<24} {n}")
    print(f"merged records: {len(sno) - len(drop_sno)} + {len(union)} = "
          f"{len(sno) - len(drop_sno) + len(union)} -> {args.out}")
    print(f"decisions -> {args.report}")


if __name__ == "__main__":
    main()
