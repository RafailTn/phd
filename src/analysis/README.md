# Cross-referencing the union against another catalogue

Two standalone scripts, independent of the `loci_extraction/` pipeline. Neither
sources `config.sh`; every path is an argument.

| script | question it answers |
|---|---|
| `snorna_overlap.py` | do the two FASTAs contain the same sequences? |
| `snorna_locate.py` | do they describe the same genomic loci? |

Run both. Sequence agreement and coordinate agreement are separate claims, and
each caught an error in the other during this analysis.

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
