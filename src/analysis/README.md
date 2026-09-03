# Cross-referencing the union against another catalogue

Three standalone scripts, independent of the `loci_extraction/` pipeline. They
do not source `config.sh`, but `paths.py` reproduces its lookup rules, so every
path has a working default and every path can still be overridden by flag --
see [Finding inputs](#finding-inputs).

| script | question it answers |
|---|---|
| `snorna_overlap.py` | do the two FASTAs contain the same sequences? |
| `snorna_locate.py` | do they describe the same genomic loci? |
| `collapse_duplicates.py` | merge them, resolving the loci both files hold |

Run the first two before the third. Sequence agreement and coordinate agreement
are separate claims, and each caught an error in the other during this
analysis.

## `snorna_overlap.py`

```
python3 scripts/analysis/snorna_overlap.py \
    --a AluACA_union_nr.fasta --b snoRNA.txt.fa \
    --align --out overlap_report.tsv
```

Four passes: exact string equality, containment (one sequence a substring of
the other), reverse-complement equality, and Smith-Waterman. Writes a per-pair
TSV; prints counts and an identity histogram.

Against `snoRNA.txt.fa` (1951 records) this gives **6 identical, 9 contained,
15 same-locus, 0 reverse-complement**.

**Report exact and containment separately.** 6 vs 15 is not a discrepancy - it
is interval-boundary convention. `hsa-novel-ACA-418` and `SNODB2065` are 113
and 114 nt at the same locus: not equal as strings, 100% identical by
alignment. Whole-string equality understates agreement between catalogues that
call the same RNA with different extents.

**There is no natural identity threshold below ~99%.** Alu-derived sequences
sit at 80-90% identity to each other by ancestry, so the histogram is a
continuum with a hump at 85-89%, not two populations. Any cut in that range is
a judgement call and should be reported as one.

**Never score similarity with `difflib.SequenceMatcher`.** Its `autojunk`
heuristic treats any element occurring in >1% of a sequence of 200+ elements as
junk; on a 4-letter alphabet that is every base, so sequences >=200 nt silently
score 0.0 and vanish from the ranking. This produced a bimodal distribution
that does not exist, and a "clean break at 80%" that was pure artifact.

`--align` needs EMBOSS `water`, which is not in `deps/pixi.toml`:

```
pixi add emboss          # or:
pixi exec -c conda-forge -c bioconda --spec emboss -- python3 ...
```

The align pass is k-mer prefiltered (`--kmer`, default 20) - aligning all 1951
records is one `water` call each and takes hours. 39 of 1951 survive the
filter. Safe above ~90% identity, lossy in the 75-85% band where hits are Alu
background anyway.

## `snorna_locate.py`

```
python3 scripts/analysis/snorna_locate.py \
    --fasta snoRNA.txt.fa --prefix SNODB \
    --genome ~/Downloads/hg38/GRCh38.primary_assembly.genome.fa \
    --bed AluACA_union_nr.bed --out work_map/snorna_placed.bed
```

Places coordinate-less sequences by exact genome search, then intersects with a
BED. Streams the assembly one record at a time (~7-10 min, memory bounded by
the largest chromosome). Same exact-matching rationale as
`loci_extraction/06_genome_search.py`.

Result: **25/25 placed uniquely, 15 overlap the union same-strand** - the same
15 the sequence comparison calls same-locus, established independently.

**Why this script exists.** `snoRNA.txt.fa` has no coordinates, only names like
`SNODB2065.id2065`. Reading a coordinate out of `snoDB_All_V2.0.tsv` by
matching that id to the `snoDB2065` row is wrong - all 25 land on a *different
chromosome* from the row of the same number; the numbering collides with
snoDB's id range by coincidence. That shortcut produced a confident and
entirely false "paralogous copies elsewhere, zero coordinate overlap"
conclusion. The tell was that every FASTA length disagreed with its assumed
coordinate span. Do not infer coordinates from an identifier that merely looks
familiar; find the sequence in the genome.

## collapse_duplicates.py

Merges `snoRNA.txt.fa` and `AluACA_union_nr.fasta` into one non-redundant
catalogue, resolving the 15 loci that both files contain.

```bash
python3 scripts/analysis/collapse_duplicates.py
```

Every path has a default -- see [Finding inputs](#finding-inputs) -- so that
bare form is usually enough. Override any of it explicitly:

```bash
python3 scripts/analysis/collapse_duplicates.py \
  --sno data/snoRNA.txt.fa --union data/AluACA_union_nr.fasta \
  --out AluACA_snoRNA_merged_nr.fasta \
  --report AluACA_snoRNA_collapse_report.tsv \
  --bed data/AluACA_union_nr.bed --bed-out AluACA_union_nr.collapsed.bed
```

A missing `--sno` or `--union` is fatal and names the three directories that
were searched. A missing `--bed` is not: the merge still runs and the
coordinate correction is skipped with a note.

All 15 pairs are exact substring relationships -- there are no internal
substitutions -- so containment finds them and no alignment is needed. Which
member survives, in order: identical -> the union record; exactly one ending
canonically (`ACA` box 3 nt from the 3' end) -> that one; otherwise the
difference is 5'-only, and a 5' difference of `--offset-tol` nt or less
(default 1) is treated as a coordinate off-by-one in favour of the union,
while a larger one keeps whichever sequence carries more of the 5' hairpin.

Expected output:

```
corrected 1 BED interval(s) -> AluACA_union_nr.collapsed.bed
pairs collapsed: 15  (union kept 14, snoDB sequence kept 1)
  identical                6
  union-ACA-canonical      4
  5prime-offset            3
  longer-5prime-hairpin    2
merged records: 1936 + 765 = 2701 -> AluACA_snoRNA_merged_nr.fasta
```

The surviving record always keeps the union header so the merged file has one
naming scheme. The single case where the snoDB *sequence* wins is
`AluACA163` / `SNODB2089`: `AluACA163` is a Jady et al. deposit and those are
3'-half only by design, so `SNODB2089` supplies 36 nt of the omitted 5'
hairpin. That extension is genomic (verified against hg38), so `--bed-out`
writes the corrected interval, chr19:23751717-23751832(-). Re-extracting the
corrected BED reproduces all 765 union sequences in the merged FASTA exactly.

One duplicate sequence survives on purpose: `hsa-novel-ACA-462` and
`hsa-novel-ACA-463` are identical over 262 nt but sit ~3 kb apart on chr19
(50101890 and 50104953, both `-`). They are two real loci, not a cataloguing
artefact, so collapsing them would lose one.


## Finding inputs

All three scripts share `paths.py`, which mirrors `loci_extraction/config.sh`
so the two halves of the project agree on where things live:

| what | looked for in, in order |
| --- | --- |
| inputs | `./`, then `$PROJ/`, then `$PROJ/data/` |
| outputs | `$OUT`, or `$PROJ` when `OUT` is unset |
| tools (`bedtools`, `water`) | `--flag`, `$BEDTOOLS`/`$WATER`, `$BIN`, `$PROJ/deps/.pixi/envs/default/bin`, `PATH` |

`$PROJ` comes from the environment when set, otherwise from walking up from
`scripts/analysis/` until a directory holds `deps/`, or holds one of the
project's marker files (`napRNA_Alu_L1_ACA.csv`, `snoRNA.txt.fa`,
`AluACA_union_nr.fasta`, `Supplemental_material.pdf`) either at its root or in
its `data/`. Several markers rather than one because a checkout may have the
raw inputs but not the outputs, or the reverse; anchoring on a single file
leaves `$PROJ` pointing at `scripts/`.

So on a server with everything under `data/`, all three run bare from
anywhere:

```bash
python3 scripts/analysis/snorna_overlap.py --align
python3 scripts/analysis/collapse_duplicates.py
```

**The one exception is the genome.** `snorna_locate.py --genome` defaults to
`$HG38_FA`, then to a lookup for `GRCh38.primary_assembly.genome.fa` in the
three directories above. A 3 GB assembly normally lives outside the project,
so that lookup usually fails and you must either export `HG38_FA` (as
`config.sh` already does) or pass `--genome`:

```bash
HG38_FA=/path/to/GRCh38.primary_assembly.genome.fa \
  python3 scripts/analysis/snorna_locate.py
```

A missing input is fatal and names the three directories searched. A missing
`--bed` in `collapse_duplicates.py` is not: the merge runs and the coordinate
correction is skipped with a note.

Note that the current directory is consulted first, so running from an
unrelated directory that happens to contain a same-named file will silently
pick that one up.
