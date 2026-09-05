# AluACA-guided chimeric reads in DKC1 chimeric eCLIP (hg38)

Runs the [VanNostrandLab/snoRNA-chimeric](https://github.com/VanNostrandLab/snoRNA-chimeric)
pipeline against the merged AluACA + snoRNA catalogue on hg38, to find chimeric reads
whose guide arm is an **AluACA**, and to say what the other arm is — rRNA, snRNA, tRNA,
or a genomic locus (which is where mRNA targets show up).

```bash
bash src/chimeric/run_all.sh          # everything: fetch, build, run, report
```

That is the whole procedure. It downloads the two runs, downloads and indexes the
references, runs the pipeline on the IP and its input control, annotates, compares
against the published output and regenerates `results/chimeric/RESULTS.md`. Each step is
skipped when its outputs are already there, so re-running is cheap, and the first run
takes roughly a day, most of it the STAR index.

`run_all.sh` takes arm names, sample accessions, or neither:

```bash
bash src/chimeric/run_all.sh arm1 arm3            # two arms of the matrix below
bash src/chimeric/run_all.sh --only report        # just regenerate RESULTS.md
bash src/chimeric/run_all.sh --species hg19 --source merged SRR30692552
```

The steps are `fastq`, `refs`, `index`, `run`, `annotate`, `compare`, `report`, selectable
with `--only a,b` or `--skip a,b`.

`report` runs only for the merged-catalogue arms. The report splits every table by
AluACA vs snoRNA guide, and the plain `snoRNA.txt.fa` catalogue holds no AluACA records,
so there would be nothing to compare against; those arms are scored by
`compare_to_published.py` instead.

There is one committed report per genome build, so the two can be read side by side:

| arm | writes |
|---|---|
| `arm0_hg38_merged` | `results/chimeric/RESULTS.md` |
| `arm3_hg19_merged` | `results/chimeric/RESULTS.hg19.md` |
| any other arm | `results/chimeric/<arm>/RESULTS.md` |

Both committed files are tracked; everything else stays out of git. Writing per-arm
means running several arms never leaves one report holding whichever finished last.

`make_report.py` follows `--gtag` throughout — the title, the genomic column of the
crosstab, the biotype and mRNA sections and the published cross-check all name the build
being reported. On an hg19 arm the cross-check additionally drops the genome-build term
from its list of expected differences, since that arm is on the published build: what is
left is the catalogue and the masking index.

Each step is also a script you can run on its own — `fetch_fastq.sh`, `fetch_refs.sh`,
`build_indices.sh`, `run_chimeras.sh`, `annotate_chimeras.py`,
`compare_to_published.py`, `make_report.py`.

`make_report.py` recomputes every number in the report from the annotated tables --
including the figures quoted in the interpretation -- so rerun it after any change rather
than editing the report by hand.

### Configuration

Every path is settable three ways, in increasing precedence: a built-in default, an
environment variable, then a command-line flag. `config.sh` holds the contract and is
sourced by every step, so a flag works on the driver and on the individual scripts alike:

```bash
bash src/chimeric/run_all.sh      --species hg19 --cpus 8
bash src/chimeric/run_chimeras.sh --source-fasta /ref/custom.fa SRR30692552
SPECIES=hg19 bash src/chimeric/build_indices.sh          # the env form still works
```

The ones worth knowing:

| flag | default | meaning |
|---|---|---|
| `--species` | `hg38` | genome tag; also the suffix of the pipeline's genomic outputs |
| `--source` | `merged` | guide catalogue: `merged`, `plain`, or a path |
| `--arm` | derived | sets species, catalogue and index density together |
| `--sparse-d` | `1` (dense) | `--sparse` builds a half-size index for a RAM-limited machine |
| `--cpus` | `nproc` | threads |
| `--proj` `--ref` `--data` `--work` `--out` | derived from the repo | relocate any part of the layout |

`bash src/chimeric/run_all.sh --help` lists all of them, including every reference and
input path. Nothing is hardcoded to one machine: the project root is found by walking up
from `config.sh`, so a checkout runs wherever it is put.

Results always land in `results/chimeric/<arm>/<SRR>/`. The default run is
`arm0_hg38_merged`.

### Running on another machine

The arms want a dense STAR index, which needs about 32 GB to build and 29 GB to
align — more than this 31 GB workstation has. `make_bundle.sh` packs a checkout small
enough to copy (~4 MB: the code, both guide catalogues, the target RNAs, the adapters and
the RepBase consensus; the genome and the FASTQs are re-fetched at the other end):

```bash
bash src/chimeric/make_bundle.sh
scp work/chimeric/chimeric-bundle.tar.gz server:
# then, there:
mkdir phd && tar -xzf chimeric-bundle.tar.gz -C phd && cd phd
pixi install --manifest-path deps/pixi.toml
bash src/chimeric/run_all.sh arm1 arm3
```

Budget >= 64 GB RAM and ~70-150 GB of disk. The bundle's file list is derived from
`config.sh` rather than restated, so an input added to the pipeline cannot be left out
of it.

## Reproducing the published run

The published output for this sample (`data/DKC1_IP.snoRNA.hg19.chimeras.csv`) holds
45,810 genomic chimeras against our 38,149 on hg38, with 28,085 read names shared.
Three things differ at once, so they are separated with a 2x2:

|          | plain `snoRNA.txt.fa` (1,951) | merged AluACA+snoRNA (2,701) |
|----------|-------------------------------|------------------------------|
| **hg38** | arm2 — catalogue effect       | arm0 — the headline run      |
| **hg19** | arm1 — the reproduction       | arm3 — build stability       |

`compare_to_published.py` scores an arm against the publication by read name (names carry
the Illumina identifier and UMI, so they are genome-independent) and traces every
published chimera the arm missed through that arm's own kept intermediates to say which
stage dropped it:

```bash
python3 src/chimeric/compare_to_published.py \
    --outdir results/chimeric/arm0_hg38_merged/SRR30692552 \
    --uid SRR30692552 --gtag hg38 --published data/DKC1_IP.snoRNA.hg19.chimeras.csv
```

For arm0 that puts **84.7 %** of the 17,725 missing reads at the genome-mapping step —
they survive masking, get a guide hit, and pass the back-map filter, then fail to yield a
chimera. Masking and catalogue competition account for the other 15 %. That is why the
build matters enough to test.

These figures come from the dense run and are regenerated with it; the committed
`*.vs_published.txt` files are the source, so re-read them rather than this paragraph
after any rerun.

Every arm needs a dense STAR index, which does not fit in this machine's 31 GB; see
[Running on another machine](#running-on-another-machine) for the bundle that runs them
elsewhere.

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
| `ref/chimeric/hg38_star_index` | GRCh38 primary + GENCODE v47, `sjdbOverhang 139`, dense |
| `ref/chimeric/guide_loci.hg38.bed` | 2533 guide loci from `snoDB_with_AluACA_union.tsv` |
| `ref/chimeric/rmsk.hg38.bed` | UCSC RepeatMasker, for the Alu-to-Alu flag |

The committed results are built on **dense** STAR indices (`--genomeSAsparseD 1`), which
is now the default and what every named arm specifies. A dense GRCh38 + sjdb index
needs ~29.4 GB resident (SA 24.9 + Genome 3.1 + SAindex 1.5) to align and more
than that to build, so it does not fit in this workstation's 31 GB; the arms behind the
current reports were run on a larger machine (see
[Running on another machine](#running-on-another-machine)).

Sparsity remains available as a fallback -- `--sparse` halves the suffix array to 12.5 GB,
giving ~17 GB resident -- and it is what the earlier results on this workstation used. A
sparse run is named `<species>_<source>_sparse` rather than an arm name, so it cannot be
mistaken for, or overwrite, a dense result. `--limitGenomeGenerateRAM` only chunks the SA
sort and cannot get peak usage below the finished SA, so it is not a substitute.

**What the sparse-to-dense move actually changed.** The same arm (hg38, merged catalogue)
was run both ways, so the density is the only variable. Scored against the published run
by read name:

| | sparse (D=2) | dense (D=1) |
|---|---|---|
| genomic chimeras called | 39,061 | 38,149 |
| shared with published | 25,005 | 28,085 |
| …as a share of this arm | 64.0% | **73.6%** |
| published-only (missed) | 20,805 | 17,725 |
| this-arm-only | 14,056 | **10,064** |

Dense calls slightly *fewer* chimeras (-2.3%) but agrees with the publication far better:
+3,080 shared reads, and this-arm-only calls fall by 28%. Most of what sparsity added was
not signal. The mechanism is that the genome step runs `--outFilterMultimapNmax 1`, so a
read is dropped unless it maps to exactly one locus, and a seed search that finds one more
or one fewer near-equal placement flips the read between kept and dropped -- an effect
concentrated in short arms.

**The enrichment ratios were unaffected, as predicted.** Both samples used the same index,
so the bias applied equally and cancelled in the IP-vs-input comparison:

| guide class | sparse | dense |
|---|---|---|
| AluACA | 0.68x (0.65 - 0.71) | 0.68x (0.65 - 0.71) |
| snoRNA | 12.75x | 12.99x |
| all | 2.60x | 2.61x |

So absolute counts and individual gene assignments moved by 1-2%, while every conclusion
in the report -- all of which rest on ratios -- is unchanged. Quote counts from the dense
run; treat per-gene assignments near the detection floor as provisional either way.

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
