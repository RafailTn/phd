# AluACA-guided chimeric reads in DKC1 chimeric eCLIP (hg38)

Runs the [VanNostrandLab/snoRNA-chimeric](https://github.com/VanNostrandLab/snoRNA-chimeric)
pipeline against the merged AluACA + snoRNA catalogue on hg38, to find chimeric reads
whose guide arm is an **AluACA**, and to say what the other arm is — rRNA, snRNA, tRNA,
or a genomic locus (which is where mRNA targets show up).

```bash
bash src/chimeric/build_indices.sh                 # once, ~1 h
bash src/chimeric/run_chimeras.sh SRR30692552      # the IP
bash src/chimeric/run_chimeras.sh SRR30692553      # the input control
python3 src/chimeric/annotate_chimeras.py --outdir results/chimeric/SRR30692552 \
    --uid SRR30692552 --out results/chimeric/SRR30692552.annotated.tsv \
    --source-bed ref/chimeric/guide_loci.hg38.bed --rmsk ref/chimeric/rmsk.hg38.bed
# ... and the same for SRR30692553, then:
python3 src/chimeric/make_report.py            # writes results/chimeric/RESULTS.md
```

`make_report.py` recomputes every number in `RESULTS.md` from the annotated tables, so
rerun it after any change rather than editing the report by hand.

`build_indices.sh` and `run_chimeras.sh` are configured by environment variable, so the
same two scripts drive every arm of the reproduction matrix below:

| variable | default | meaning |
|---|---|---|
| `SPECIES` | `hg38` | genome tag; also the suffix of the pipeline's genomic outputs |
| `GENOME_INDEX` | `ref/chimeric/${SPECIES}_star_index` | STAR index to use |
| `SOURCE_FASTA` | `data/AluACA_snoRNA_merged_nr.fasta` | the guide catalogue |
| `SPARSE_D` | `2` | `build_indices.sh` only; `1` builds a dense index |

## Reproducing the published run

The published output for this sample (`data/DKC1_IP.snoRNA.hg19.chimeras.csv`) holds
45,810 genomic chimeras against our 39,061 on hg38, and only 25,005 read names are
shared. Three things differ at once, so they are separated with a 2x2:

|          | plain `snoRNA.txt.fa` (1,951) | merged AluACA+snoRNA (2,701) |
|----------|-------------------------------|------------------------------|
| **hg38** | arm2 — catalogue effect       | arm0 — the headline run      |
| **hg19** | arm1 — the reproduction       | arm3 — build stability       |

`compare_to_published.py` scores an arm against the publication by read name (names carry
the Illumina identifier and UMI, so they are genome-independent) and traces every
published chimera the arm missed through that arm's own kept intermediates to say which
stage dropped it:

```bash
python3 src/chimeric/compare_to_published.py --outdir results/chimeric/SRR30692552 \
    --uid SRR30692552 --gtag hg38 --published data/DKC1_IP.snoRNA.hg19.chimeras.csv
```

For arm0 that puts **85.3 %** of the 20,805 missing reads at the genome-mapping step —
they survive masking, get a guide hit, and pass the back-map filter, then fail to yield a
chimera. Masking and catalogue competition account for the other 15 %. That is why the
build matters enough to test.

The hg19 arms need a dense STAR index, which does not fit in this machine's 31 GB; see
[`server/README.md`](server/README.md) for the bundle that runs them elsewhere.

## Which SRA run to use

`data/` holds four runs. Only two are chimeric eCLIP:

| run | GSM | reads | length | what |
|---|---|---|---|---|
| **SRR30692552** | GSM8521923 | 14.85 M | 150 | **DKC1 IP chimeric eCLIP** |
| SRR30692553 | GSM8521922 | 3.48 M | 150 | DKC1 chimeric eCLIP input |
| SRR30693619 | GSM8521973 | 19.74 M | 75 | standard DKC1 eCLIP |
| SRR30693620 | GSM8521972 | 23.44 M | 75 | standard DKC1 eCLIP |

**The NCBI SRA web search page pairs these titles with the wrong runs** — it swaps the
two groups. `efetch -db sra -rettype runinfo` has it right, and so does the direct
evidence: all 45,810 reads in `data/DKC1_IP.snoRNA.hg19.chimeras.csv` — the published
output for GSM8521923 — carry flowcell `NB552078:290:HYGLFAFXY`, which is SRR30692552's,
and their sequences run to 130 nt = 150 − 10 UMI − 10 randomer, impossible for a 75 nt run.

That published CSV is a genuine ground truth: same sample, just hg19 and the plain
1951-record `snoRNA.txt.fa` source instead of hg38 and the merged catalogue.

## The adapter file, which is not published

`sno-chimeras.py --strategy chimeric` reads
`/storage/vannostrand/software/eclip/data/se.2.round.adapters.fasta`, a hard-coded lab
path. The file is in no public repo, so it was reconstructed from the published output
and lives at `ref/chimeric/se.2.round.adapters.fasta`.

It is **24 tiles of 10-mers** sliding across `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`
(the RiL19 / TruSeq read-through), built the way
`VanNostrandLab/clip/eclip/bin/generate_adaptertrim_fasta.ipynb` builds its InvRNA
files — `for i in range(len(seq) - word_len)`, so **the last tile is dropped** — but
with `word_len = 10` rather than that notebook's 15.

With the pipeline's own `--times 3 -e 0.1 -q 6 -O 1`, a 10 nt 5' UMI and `cutadapt
-u -10`, this reproduces **4983 / 5000 (99.66 %)** of the published trimmed sequences
byte-for-byte.

<details><summary>How it was identified, if it ever needs redoing</summary>

Every published `sequence` is an exact prefix of its raw read after removing 10 nt, and
the *minimum* removal is exactly 13 nt = 10 (`-u -10`) + 3. That 3 is the tell: the tile
set's first bases cover all of A/C/G/T, so with `-O 1` every read loses one base per
round, and `--times 3` loses three. Sweeping `--times` (1/2/3/4) × `-q` (6/15/20) gives
a sharp optimum at `times=3 q=6` (221/300 vs ≤33 elsewhere); sweeping tile length then
peaks hard at 10 (299/300 vs 221 at 15).

Plausible-looking wrong answers: the 15-mer RiL19 tiles from `YeoLab/eclip` (221/300),
the concatenated InvRNA1–8 tiles (6/300 — their leading `NN` plus
`--match-read-wildcards` matches everything), anything including the TruSeq index region
(1/300). Cutadapt version is *not* the variable — 3.2 and 5.2 score within one read of
each other, so no version pin is needed.
</details>

## Why the source catalogue keeps the snoRNAs in

`--source_rna_fasta` is `data/AluACA_snoRNA_merged_nr.fasta`: the 765-locus AluACA union
merged with the 1936 records of `snoRNA.txt.fa`. Running AluACAs alone would be a
mistake — bowtie2 scores each read against the whole catalogue at once, so keeping the
real snoRNAs in means an AluACA call is an AluACA *beating every snoRNA*, not an AluACA
being the only thing on offer.

AluACA records are exactly those whose FASTA id ends `.id3xxx` (3001–3765, numbered above
`snoRNA.txt.fa`'s highest id so the two sets cannot collide).

On a 321k-read pilot this separates cleanly: of reads with a best-scoring guide,
**5.9 % best-hit an AluACA, 94.0 % a snoRNA, and only 0.1 % are ambiguous between the
two**. Adding the 765 AluACAs costs only 3.4 % more aligned reads and takes mean
alignments per read from 7.4 to 9.4 — the Alu redundancy does not blow up `bowtie2 -a`.

## Changes to the upstream code

`sno-chimeras.py` and `scripts/` are vendored verbatim except for these, all marked in
`git log`:

1. **`convert_sam_to_bam` crashed on `--source_rna_bed`** — `str.replace()` was called
   with cmder keyword arguments, and it then indexed an undefined `true_bam`. Fixed, so
   the guide-locus filter is usable.
2. **`individual_chimeras` shelled out to `clipper`** for the genome tag unconditionally.
   clipper is not installed, has no hg38 model, and the authors had already commented out
   their own genomic peak-calling task; peak calling is now gated on
   `--clipper_rna_species`.
3. **Hard-coded lab paths** for the adapter FASTAs now resolve under `$CHIMERIC_REF_DIR`,
   and the emitted run-log no longer sources a cluster-only venv.
4. `logger.error` in the strategy else-branch referenced an undefined `strategy`.

Two console scripts the pipeline calls come from the lab's unpublished `iToolBox` and are
reimplemented in `scripts/`: `fastq_to_fasta` (note it takes two bare positional
arguments, *not* the FASTX-Toolkit `-i/-o` interface, so the FASTX binary cannot be
substituted) and `sno_chimeras_summary` (command-line contract only, not the original
report layout).

The pipeline runs with `--keep`, because its cleanup step deletes the per-target
`.chimeras.csv` files that are the actual result.

## References built here

| path | what |
|---|---|
| `ref/chimeric/se.2.round.adapters.fasta` | reconstructed adapter tiles (above) |
| `ref/chimeric/repbase/human_repbase.fa` | 1144 consensus sequences, RepBase `humrep.ref` + `humsub.ref` |
| `ref/chimeric/repbase_star_index` | STAR index of the above |
| `ref/chimeric/hg38_star_index` | GRCh38 primary + GENCODE v47, `sjdbOverhang 139` |
| `ref/chimeric/guide_loci.hg38.bed` | 2533 guide loci from `snoDB_with_AluACA_union.tsv` |
| `ref/chimeric/rmsk.hg38.bed` | UCSC RepeatMasker, for the Alu-to-Alu flag |

`hg38_star_index` is built with `--genomeSAsparseD 2`. GRCh38 + sjdb does not otherwise
fit in this machine's 31 GB: a dense index needs ~29.4 GB resident (SA 24.9 + Genome 3.1 +
SAindex 1.5) against ~26 GB available, which is why a default `genomeGenerate` gets
OOM-killed here. Sparsity halves the suffix array to 12.5 GB, giving ~17 GB resident.
`--limitGenomeGenerateRAM` only chunks the SA sort and cannot get peak usage below the
finished SA, so it is not a substitute.

**Sparsity is not results-neutral here, despite the manual describing it as a speed
tradeoff.** Measured on a controlled A/B -- dense and sparse indices of chr1, identical
in every other parameter, queried with this pipeline's own target arms under its own
genome-step settings:

| | dense (D=1) | sparse (D=2) |
|---|---|---|
| uniquely mapped | 20,600 | 19,134 |
| mapped by both | 15,456 | |
| …same chrom+pos+CIGAR | 14,482 (93.7%) | |
| …placed differently | 974 (6.3%) | |
| lost by sparse / gained by sparse | 5,144 | 3,678 |

The instability concentrates in short arms -- 13.7% disagreement at 16-20 nt, 9.6% at
21-25 nt, ~2% above 40 nt. The mechanism is not that sparse alignments are wrong but that
this step runs `--outFilterMultimapNmax 1`, so a read is discarded unless it maps to
exactly one locus; a seed search that finds one more or one fewer near-equal placement
flips the read between kept and dropped.

Two things follow for how the results should be read:

- **IP-vs-input folds are robust.** Both samples use the same index and their target-arm
  length distributions are nearly identical (median 43 nt each; 22.2% vs 21.0% in the
  unstable 16-25 nt band), so the bias applies equally and cancels in the ratio.
- **Absolute counts and individual gene assignments are not.** Expect ~7% uncertainty in
  totals, and more for snoRNA-guided chimeras than AluACA-guided ones, since snoRNA arms
  are shorter (median 30 nt, 43.6% in the unstable band) than AluACA arms (median 46 nt,
  23.1%).

The chr1-only design inflates these numbers -- a read whose true locus is on another
chromosome is forced onto a chr1 paralog -- so treat them as an upper bound rather than an
estimate of the effect on the real hg38 run. Rebuilding dense on a >32 GB machine would
remove the caveat entirely; that is `arm0d` in
[`server/README.md`](server/README.md), which measures the same thing on real chimera
calls genome-wide and supersedes this test once it lands.

<details><summary>Redoing the A/B, if it is ever needed</summary>

`work/chimeric/sparsetest/` held the indices, the query and the two alignments. It was
deleted: the indices and `chr1.fa` are regenerable, and the query was a byte-for-byte copy
of a file the run directory still holds. The whole thing rebuilds in about five minutes.

The query is the 156,616 target arms that reached the genome-mapping step -- i.e. what
survived the back-map filter and did not hit rRNA, snRNA or tRNA. **It only exists inside
a completed run**, so this cannot be redone if `results/chimeric/SRR30692552/` is cleared.

```bash
cd /home/rafail/Downloads/phd
export PATH=$PWD/deps/.pixi/envs/default/bin:$PATH
S=work/chimeric/sparsetest; mkdir -p $S

cp results/chimeric/SRR30692552/SRR30692552.snoRNA.RNA.unmap.fasta $S/query.fasta
samtools faidx /home/rafail/Downloads/hg38/GRCh38.primary_assembly.genome.fa chr1 > $S/chr1.fa

for D in 1 2; do
  mkdir -p $S/idx_D$D $S/aln_D$D
  STAR --runMode genomeGenerate --runThreadN 20 --genomeDir $S/idx_D$D \
       --genomeFastaFiles $S/chr1.fa --genomeSAindexNbases 13 \
       --genomeSAsparseD $D --outFileNamePrefix $S/idx_D$D/
  # These are the pipeline's own genome-step settings, which is what makes the
  # comparison meaningful -- in particular --outFilterMultimapNmax 1, the filter
  # the sparsity interacts with.
  STAR --alignEndsType EndToEnd --genomeDir $S/idx_D$D --genomeLoad NoSharedMemory \
       --outFileNamePrefix $S/aln_D$D/ \
       --outFilterMatchNminOverLread 0.66 --outFilterMultimapNmax 1 \
       --outFilterMultimapScoreRange 1 --outFilterScoreMin 10 \
       --outFilterScoreMinOverLread 0.66 --outFilterType BySJout \
       --outReadsUnmapped Fastx --outSAMattributes Standard --outSAMmode Full \
       --outSAMtype SAM --outSAMunmapped None --outStd Log \
       --readFilesIn $S/query.fasta --runMode alignReads --runThreadN 20
done
```

Then compare the two `Aligned.out.sam`, ignoring secondary alignments (flag `0x100`), on
`(reference, position, CIGAR)` per read name.
</details>

The repeat index is built from RepBase human consensus only. The original run used the
lab's `homo_sapiens_repbase_v2`, which is not public and probably also carried rRNA and a
few other abundant sequences, so the masking stage here is slightly less aggressive than
theirs.

## Does masking discard real AluACA signal?

A fair worry, since every AluACA is Alu-derived and so is much of RepBase. Both masking
stages were tested by taking the reads each one discarded and pushing them through the
rest of the pipeline to see how many would have become AluACA chimeras.

**Repeat masking removes almost nothing, and loses almost nothing.**

| | reads |
|---|---|
| discarded by repeat masking | 185,300 (1.57% of input) |
| …would have died at genome masking anyway | 165,672 |
| …uniquely lost to the repeat step | 19,628 |
| …yielding any putative chimeric target | 101 |
| …with an AluACA guide | **72** |

72 against the 13,759 AluACA chimeras called is **0.52%**, and that is an upper bound --
those 72 still had to survive back-mapping and the gap/length criteria.

The reason is `--alignEndsType EndToEnd`. Soft-clipping is forbidden, so the whole read
must align, and `--outFilterMismatchNmax 10` caps mismatches at 10. In the called AluACA
chimeras the non-guide arm has a median of 54 nt and **100%** exceed 10 nt, so no chimera
can be forced onto an Alu consensus within the mismatch budget. A read that is Alu along
its entire length does get masked -- which is the point, since that is a non-chimeric
AluACA read. The thresholds are also mutually consistent: an arm short enough to slip
under the 10-mismatch budget is already below the >=16 nt target minimum and could never
have become a chimera.

**Genome masking removes a great deal, and it is right to.** It discards 7,226,001 reads
(62%), among which 29,872 would have been called AluACA chimeras -- 2.2x the number
actually reported. But classifying how STAR aligned them settles what they are:

| how STAR aligned the discarded read | n | share |
|---|---|---|
| contiguous (CIGAR all M) | 29,321 | 98.2% |
| spliced across an annotated GENCODE junction | 101 | 0.3% |
| spliced across a **novel** junction | 450 | 1.5% |

98.2% align as one unbroken stretch of genome. They are single transcripts, not chimeras:
the chimeric appearance comes from the source mapper calling an Alu-like segment as an
AluACA guide and the sequence next to it as a "target". Removing them is the whole purpose
of the stage -- without it the AluACA set would be more than tripled with artefacts.

Only the 450 novel-junction reads (median gap 1,505 bp) are genuinely ambiguous: STAR
joined two blocks across a gap that is not an annotated splice site, which is what a real
*cis* interaction would also look like, since STAR permits introns up to ~590 kb by
default. That is **3.3%** of the called AluACA set, and an upper bound, since unannotated
splicing accounts for an unknown share of it.

**Combined upper bound on masking loss: ~522 reads, under 4% of the AluACA chimeras.**
The one loss mode neither test can rule out is an AluACA paired with an *exonized Alu*,
where both arms are Alu-like and the read can be spliced-aligned within the Alu consensus.
That case is unrecoverable here -- and also indistinguishable from the Alu-to-Alu artefact
that is already 37.9% of what survives.

## Reading the output

Per-target `<uid>.snoRNA.<target>.chimeras.csv` from the pipeline, pooled and annotated
by `annotate_chimeras.py` into one TSV with:

- `guide_class` — `AluACA` / `snoRNA` / `ambiguous`
- `target_class` — `rRNA` / `snRNA` / `tRNA` / `hg38`
- `gene_name`, `gene_type`, `feature` — GENCODE annotation of the genomic arm.
  **`gene_type == protein_coding` with `feature == exonic` is an AluACA–mRNA chimera.**
- `target_in_source_locus` — the genomic arm lands on a guide locus, same strand, ≥25 %
  overlap. Usually one contiguous transcript rather than a real chimera.
- `target_in_repeat` — the genomic arm lands in an annotated repeat. Every AluACA is
  Alu-derived and Alu is the largest repeat family in the genome, so **Alu-to-Alu is the
  dominant false-positive mode for this particular question** and this column is the one
  to look at before believing any AluACA target.

Both flags are columns, not filters, so the cost of applying either stays visible.
