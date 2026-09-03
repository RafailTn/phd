#!/usr/bin/env python3
"""
Place the sequences of a coordinate-less FASTA in hg38 by exact search, then
intersect with a BED to see whether the two sets describe the same loci.

Why this exists: snoRNA.txt.fa carries no coordinates, only names of the form
`SNODB2065.id2065`. The tempting shortcut is to read a coordinate out of
snoDB_All_V2.0.tsv by matching that id to the `snoDB2065` row. That mapping is
WRONG - all 25 SNODB entries land on a different chromosome from the snoDB row
of the same number - and using it produced a confident, entirely false "these
are paralogous copies elsewhere in the genome, zero coordinate overlap"
conclusion. The names collide with snoDB's id range by coincidence.

So: never infer coordinates from an identifier that happens to look familiar.
Find the sequence in the genome.

Exact matching is deliberate and is the same approach as
loci_extraction/06_genome_search.py: despite being Alu-derived, individual
copies carry enough substitutions that a ~100 nt exact match is usually unique.
Every hit is reported, so a sequence matching several loci stays visibly
ambiguous instead of being silently collapsed to the first.

Streams the genome one record at a time, so memory is bounded by the largest
chromosome (~250 MB), not the 3.1 GB assembly. Takes roughly 7-10 minutes.

Usage:
  python3 snorna_locate.py --fasta snoRNA.txt.fa --bed AluACA_union_nr.bed \
                           --genome ~/Downloads/hg38/GRCh38.primary_assembly.genome.fa \
                           --out work_map/snorna_placed.bed
  # restrict to a subset of records:
  python3 snorna_locate.py --fasta snoRNA.txt.fa --prefix SNODB ...

Output BED is 0-based half-open, matching AluACA_union_nr.bed.
"""
import argparse
import os
import subprocess
import shutil
import sys

from paths import find_input, find_tool, out_path, require

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path, prefix=None):
    seqs, name, buf = {}, None, []

    def flush():
        if name is not None and (prefix is None or name.startswith(prefix)):
            seqs[name] = "".join(buf).upper()

    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    flush()
    return seqs


def find_all(hay, needle):
    i = hay.find(needle)
    while i != -1:
        yield i
        i = hay.find(needle, i + 1)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", default=find_input("snoRNA.txt.fa"),
                   help="FASTA to place (no coordinates)")
    p.add_argument("--genome",
                   default=os.environ.get("HG38_FA")
                           or find_input("GRCh38.primary_assembly.genome.fa"),
                   help="genome FASTA (uncompressed); $HG38_FA is honoured")
    p.add_argument("--bed", default=find_input("AluACA_union_nr.bed"),
                   help="BED to intersect the placements against")
    p.add_argument("--out", default=out_path("work_map/snorna_placed.bed"),
                   help="output BED of placements")
    p.add_argument("--prefix", help="only place records whose name starts with this")
    p.add_argument("--bedtools", default=None)
    args = p.parse_args()

    require([("FASTA to place (--fasta)", args.fasta),
             ("genome FASTA (--genome)", args.genome)])
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    queries = read_fasta(args.fasta, args.prefix)
    if not queries:
        sys.exit("no records selected from " + args.fasta)
    print(f"placing {len(queries)} sequence(s) from {os.path.basename(args.fasta)}")

    hits = {k: [] for k in queries}

    def scan(chrom, seq):
        if chrom is None:
            return
        for k, v in queries.items():
            for pat, strand in ((v, "+"), (rc(v), "-")):
                for i in find_all(seq, pat):
                    hits[k].append((chrom, i, i + len(v), strand))

    chrom, buf = None, []
    for line in open(args.genome):
        if line[0] == ">":
            scan(chrom, "".join(buf).upper())
            chrom, buf = line[1:].split()[0], []
            print(f"\r  scanning {chrom} ...      ", end="", file=sys.stderr)
        else:
            buf.append(line.rstrip("\n"))
    scan(chrom, "".join(buf).upper())
    print("", file=sys.stderr)

    placed = sum(1 for v in hits.values() if v)
    multi = [k for k, v in hits.items() if len(v) > 1]
    missing = [k for k, v in hits.items() if not v]
    print(f"  exactly placed:      {placed}/{len(queries)}")
    print(f"  multiple exact hits: {len(multi)}" + (f"  {multi}" if multi else ""))
    print(f"  not found:           {len(missing)}" + (f"  {missing}" if missing else ""))
    if missing:
        print("    (a sequence absent from the assembly is not genomic as given - "
              "it may be padded, trimmed, or from another build)")

    rows = sorted((c, s, e, k, st) for k, v in hits.items() for c, s, e, st in v)
    with open(args.out, "w") as fh:
        for c, s, e, k, st in rows:
            fh.write(f"{c}\t{s}\t{e}\t{k}\t0\t{st}\n")
    print(f"  wrote {len(rows)} placement(s) to {args.out}")

    if not args.bed:
        return
    bedtools = find_tool("bedtools", args.bedtools)
    if not bedtools:
        sys.exit("bedtools not found; pass --bedtools PATH to run the intersect")

    for label, flag in (("same-strand", ["-s"]), ("either-strand", [])):
        r = subprocess.run([bedtools, "intersect", "-a", args.out, "-b", args.bed,
                            "-wa", "-wb"] + flag,
                           capture_output=True, text=True, check=True)
        lines = [l for l in r.stdout.splitlines() if l.strip()]
        print(f"\n{label} overlaps with {os.path.basename(args.bed)}: "
              f"{len(lines)}/{len(rows)}")
        if flag:
            for l in lines:
                f = l.split("\t")
                print(f"  {f[3]:20s} {f[0]}:{f[1]}-{f[2]}({f[5]})"
                      f"  <->  {f[9]} {f[6]}:{f[7]}-{f[8]}({f[11]})")


if __name__ == "__main__":
    main()
