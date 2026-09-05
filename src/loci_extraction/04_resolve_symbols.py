#!/usr/bin/env python3
"""
Step 04 - Map each AluACA to a current (GENCODE v47) gene symbol.

Supplemental Table 3 is from 2012, so ~66 of its symbols are deprecated
aliases (LEPRE1->P3H1, MLL3->KMT2C, BAT3->BAG6, ...). Resolution order:

  1. MANUAL_DESC  - entries whose "symbol" field is actually prose, where the
                    real symbol is recoverable by reading the description.
  2. direct       - symbol already valid in GENCODE v47.
  3. NCBI Gene    - esearch [Gene Name] + esummary -> nomenclaturesymbol,
                    accepted only if the result is itself in GENCODE.
  4. unresolved   - cDNA/EST accessions and truncated symbols with no modern
                    equivalent; these fall through to the genome-wide search
                    in step 06.

"No apparent host gene" rows (12 of them) are marked NA by design - the paper
asserts these have no host, and step 06 confirms most land intergenic.

Writes: work_map/id_gene.tsv  (aluaca_id, gene, source)
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

WORK = os.environ["WORK"]
# config.sh exports both; hardcoding them here meant --acc-from and friends
# reached every step except this one.
EUTILS = os.environ.get("EUTILS", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils")
TOOL = os.environ.get("NCBI_TOOL", "claude_code")

# Symbol buried in prose rather than in the symbol field.
MANUAL_DESC = {
    "AluACA12":  "METTL10",    # "methyltransferase like 10"
    "AluACA13":  "C19orf70",   # "...open reading frame 70 (C19orf70)"
    "AluACA14":  "HSCB",       # "HscB iron-sulfur cluster co-chaperone homolog"
    "AluACA170": "TBCE",       # "TBCE< beta-tubulin cofactor E"  (stray '<')
    "AluACA219": "LOC129293",  # "hypothetical protein LOC129293 precursor"
    "AluACA274": "DDX50",      # "DDx50" (lowercase x)
}
# Recovered from the description after NCBI returned nothing.
MANUAL_POST = {
    "AluACA68":  "RPRD2",      # "KPRD2" is a typo for RPRD2
    "AluACA280": "SESTD1",     # "SEC14 and spectrin domains 1"
    "AluACA323": "AGTPBP1",    # "KIAA1035, ... ATP/GTP binding protein 1"
}


def get(url, tries=4):
    for _ in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return r.read().decode()
        except Exception:
            time.sleep(1.5)
    return ""


def ncbi_current_symbol(sym):
    """Return the current official symbol for a legacy alias, or ''."""
    q = urllib.parse.quote(f"{sym}[Gene Name] AND Homo sapiens[Organism]")
    r = get(f"{EUTILS}/esearch.fcgi?db=gene&term={q}&retmode=json&retmax=5&tool={TOOL}")
    time.sleep(0.4)  # NCBI: 3 req/s without an API key
    try:
        ids = json.loads(r)["esearchresult"].get("idlist", [])
    except Exception:
        ids = []
    if not ids:
        return ""
    s = get(f"{EUTILS}/esummary.fcgi?db=gene&id={ids[0]}&retmode=json&tool={TOOL}")
    time.sleep(0.4)
    try:
        d = json.loads(s)["result"][ids[0]]
        return d.get("nomenclaturesymbol") or d.get("name") or ""
    except Exception:
        return ""


def main():
    gencode = {l.strip() for l in open(os.path.join(WORK, "gene_names.txt"))}
    rows = [l.rstrip("\n") for l in open(os.path.join(WORK, "table3.txt"))]

    staged, pending = [], []
    for line in rows:
        parts = [p.strip() for p in line.split(",")]
        aid = parts[0]
        sym = parts[1] if len(parts) > 1 else ""
        if aid in MANUAL_DESC:
            sym = MANUAL_DESC[aid]
        if re.search(r"no apparent host gene", line, re.I):
            staged.append((aid, "NA", "none"))
        elif sym in gencode:
            staged.append((aid, sym, "gencode_direct"))
        else:
            pending.append((aid, sym))
            staged.append((aid, sym, "PENDING"))

    sys.stderr.write(
        f"direct={sum(1 for s in staged if s[2]=='gencode_direct')} "
        f"none={sum(1 for s in staged if s[2]=='none')} pending={len(pending)}\n")

    cache = {}
    for aid, sym in pending:
        if sym not in cache:
            cache[sym] = ncbi_current_symbol(sym)
            sys.stderr.write(f"  {sym} -> {cache[sym] or 'UNRESOLVED'}\n")

    final = []
    for aid, sym, src in staged:
        if src != "PENDING":
            final.append((aid, sym, src))
            continue
        new = cache.get(sym, "")
        if new and new in gencode:
            final.append((aid, new, f"ncbi_alias:{sym}"))
        elif aid in MANUAL_POST:
            final.append((aid, MANUAL_POST[aid], "manual_desc"))
        elif new:
            final.append((aid, new, f"ncbi_notin_gencode:{sym}"))
        else:
            final.append((aid, sym, "UNRESOLVED"))

    # MANUAL_POST may also apply to rows NCBI resolved to something wrong.
    final = [(a, MANUAL_POST[a], "manual_desc") if a in MANUAL_POST else (a, g, s)
             for a, g, s in final]

    out = os.path.join(WORK, "id_gene.tsv")
    with open(out, "w") as f:
        f.write("aluaca_id\tgene\tsource\n")
        for a, g, s in final:
            f.write(f"{a}\t{g}\t{s}\n")

    counts = {}
    for _, _, s in final:
        k = s.split(":")[0]
        counts[k] = counts.get(k, 0) + 1
    for k, v in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {k:22s} {v}")
    print(f"[04] wrote {out}")


if __name__ == "__main__":
    main()
