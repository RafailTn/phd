#!/usr/bin/env python3
"""
Helper for step 05 - exact search of each AluACA sequence inside the genomic
span of its annotated host gene.

Searches both the sequence and its reverse complement, so the reported strand
is derived from the match rather than assumed. Coordinates are converted back
to genome space by adding the gene span's start offset (BED, 0-based).

Writes: work_map/located.tsv
"""
import os
import sys

WORK = os.environ["WORK"]
FASTA = os.environ["FASTA"]

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    """Yield ((aluaca_id, accession), sequence). ID comes from the defline."""
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            parts = line[1:].split()
            acc = parts[0]
            aid = [p for p in parts if p.startswith("AluACA")]
            name = (aid[0] if aid else acc, acc)
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def main():
    seqs = read_fasta(FASTA)

    gene = {}
    for i, l in enumerate(open(os.path.join(WORK, "id_gene.tsv"))):
        if i == 0:
            continue
        a, g, s = l.rstrip("\n").split("\t")
        gene[a] = (g, s)

    # gene symbol -> list of (chrom, span_start, plus_strand_sequence)
    loci = {}
    for l in open(os.path.join(WORK, "gene_seqs.tsv")):
        f = l.rstrip("\n").split("\t")
        loci.setdefault(f[3], []).append((f[0], int(f[1]), f[6].upper()))

    out = open(os.path.join(WORK, "located.tsv"), "w")
    out.write("aluaca_id\taccession\tgene\tgene_src\tchrom\tstart\tend\t"
              "strand\tseq_len\tn_hits\tstatus\n")

    stat = {}
    for (aid, acc), s in sorted(seqs.items(),
                                key=lambda x: int(x[0][0].replace("AluACA", ""))):
        g, src = gene.get(aid, ("NA", "missing"))
        hits = []
        if g in loci:
            r = rc(s)
            for chrom, gs, gseq in loci[g]:
                i = gseq.find(s)
                while i >= 0:
                    hits.append((chrom, gs + i, gs + i + len(s), "+"))
                    i = gseq.find(s, i + 1)
                i = gseq.find(r)
                while i >= 0:
                    hits.append((chrom, gs + i, gs + i + len(r), "-"))
                    i = gseq.find(r, i + 1)
        if hits:
            st = "unique" if len(hits) == 1 else "multi_in_gene"
            for c, a1, b1, strand in hits:
                out.write(f"{aid}\t{acc}\t{g}\t{src}\t{c}\t{a1}\t{b1}\t"
                          f"{strand}\t{len(s)}\t{len(hits)}\t{st}\n")
        else:
            st = "no_gene" if (g == "NA" or src == "UNRESOLVED") else "not_found_in_gene"
            out.write(f"{aid}\t{acc}\t{g}\t{src}\t-\t-\t-\t-\t{len(s)}\t0\t{st}\n")
        stat[st] = stat.get(st, 0) + 1
    out.close()

    for k, v in sorted(stat.items(), key=lambda x: -x[1]):
        print(f"  {k:22s} {v}")


if __name__ == "__main__":
    main()
