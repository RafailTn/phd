#!/usr/bin/env python3
"""
Step 06 - Genome-wide exact search for the AluACAs step 05 could not place.

Covers two cases:
  no_gene            - 12 "No apparent host gene" rows + 6 whose 2012 "symbol"
                       was a cDNA/EST accession with no modern equivalent.
  not_found_in_gene  - the sequence was absent from the assigned gene's span
                       (wrong/outdated assignment, or a neighbouring locus).

Exact matching is deliberate: despite these being Alu-derived, each deposited
copy carries enough individual substitutions that an exact 76-80 nt match is
usually unique. Sequences that still hit multiple loci are genuine paralogs
(NBPF family, chrX tandem arrays) and are reported with every hit so the
ambiguity stays visible rather than being silently collapsed.

Streams the genome one record at a time to bound memory at ~1 chromosome.

Writes: work_map/genome_hits.tsv
"""
import os
import sys

WORK = os.environ["WORK"]
FASTA = os.environ["FASTA"]
HG38 = os.environ["HG38_FA"]

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def main():
    want = set()
    for i, l in enumerate(open(os.path.join(WORK, "located.tsv"))):
        if i == 0:
            continue
        f = l.rstrip("\n").split("\t")
        if f[10] in ("no_gene", "not_found_in_gene"):
            want.add(f[0])

    seqs, name, buf = {}, None, []
    for line in open(FASTA):
        if line.startswith(">"):
            if name and name in want:
                seqs[name] = "".join(buf).upper()
            p = line[1:].split()
            a = [x for x in p if x.startswith("AluACA")]
            name = a[0] if a else p[0]
            buf = []
        else:
            buf.append(line.strip())
    if name and name in want:
        seqs[name] = "".join(buf).upper()

    sys.stderr.write(f"searching {len(seqs)} sequences genome-wide\n")

    pats = ([(a, s, "+") for a, s in seqs.items()] +
            [(a, rc(s), "-") for a, s in seqs.items()])
    hits = {a: [] for a in seqs}

    def scan(chrom, seq):
        if not seq:
            return
        for aid, pat, strand in pats:
            i = seq.find(pat)
            while i >= 0:
                hits[aid].append((chrom, i, i + len(pat), strand))
                i = seq.find(pat, i + 1)

    chrom, buf = None, []
    for line in open(HG38):
        if line.startswith(">"):
            if chrom:
                scan(chrom, "".join(buf).upper())
            chrom = line[1:].split()[0]
            buf = []
            sys.stderr.write(f"  {chrom}\r")
        else:
            buf.append(line.strip())
    if chrom:
        scan(chrom, "".join(buf).upper())
    sys.stderr.write("\n")

    out = os.path.join(WORK, "genome_hits.tsv")
    with open(out, "w") as o:
        o.write("aluaca_id\tchrom\tstart\tend\tstrand\tn_hits\n")
        for a in sorted(hits, key=lambda x: int(x.replace("AluACA", ""))):
            h = hits[a]
            if not h:
                o.write(f"{a}\t-\t-\t-\t-\t0\n")
            for c, s, e, st in h:
                o.write(f"{a}\t{c}\t{s}\t{e}\t{st}\t{len(h)}\n")
    print(f"[06] wrote {out}")


if __name__ == "__main__":
    main()
