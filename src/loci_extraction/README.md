# AluACA → hg38 mapping pipeline

Places the 348 AluACA RNA genes deposited by Jády, Ketele & Kiss
(EMBL `HE855917`–`HE856264`) onto hg38, then compares them to the napRNAdb
Alu/L1 ACA loci in `napRNA_Alu_L1_ACA.csv`.

Run everything with:

```bash
bash scripts/run_all.sh
```

## Outputs (written to `--out`, default the project root)

| file | contents |
|---|---|
| `AluACA_HE855917-HE856264.fasta` | the 348 deposited sequences (AluACA1–348) |
| `AluACA_hg38_coordinates.tsv` | 344 located AluACAs: ID, accession, host gene, hg38 coords, strand, method, matching CSV entry |
| `AluACA_hg38.bed` | the same as BED6 |
| `AluACA_unresolved.tsv` | the 4 that could not be placed unambiguously |
| `AluACA_union_nr.{bed,tsv,fasta}` | the 765-locus union set (see below) |

## Steps

| script | does |
|---|---|
| `01_fetch_sequences.sh` | NCBI E-utilities fetch + integrity checks |
| `02_extract_table3.sh` | parse host genes out of the supplemental PDF |
| `03_build_gencode_genes.sh` | flatten GENCODE v47 genes to BED |
| `04_resolve_symbols.py` | map 2012-era gene symbols to current ones |
| `05_locate_in_genes.sh` + `lib_locate.py` | exact search inside each host gene |
| `06_genome_search.py` | genome-wide exact search for the leftovers |
| `07_assemble_and_match.sh` | merge, compare to the CSV, write deliverables |
| `08_build_union.sh` | non-redundant union with the napRNAdb loci |
| `09_add_repeat_family.sh` | annotate each locus with its RepeatMasker element |

`csv_overlap_check.sh` is standalone and unrelated to the mapping — it checks
whether the two `napRNA_Alu_L1_*.csv` files overlap.

## The central design decision

These sequences are 71–162 nt (median 79) fragments of **Alu consensus**. A
genome-wide aligner returns hundreds of equally good hits per sequence,
because hg38 contains over a million Alu copies — alignment alone cannot
place them.

What makes the mapping tractable is Supplemental Table 3, which names the
**host gene** for each AluACA. The pipeline therefore restricts the search
space to that gene's span and requires an *exact* match. That yields 325
unique placements. Only the 23 sequences with no usable gene assignment go to
a genome-wide search, which resolves 19 more uniquely.

Two independent signals say the approach is sound:

- **AluACA55** was assigned only a truncated `MBOAT` in the table, yet the
  genome-wide search independently placed it inside **MBOAT1**.
- **10 of the 12** "No apparent host gene" entries land **intergenic**,
  matching the paper's own claim.

## Results

**344 / 348 placed.** 325 gene-anchored, 18 genome-unique, 1 disambiguated.

**Match to `napRNA_Alu_L1_ACA.csv`: 122 of 344** (and 121 of the 543 CSV loci
carry an AluACA). The matches are exact — they share interval boundaries with
the CSV entries, which also confirms the CSV is hg38. Since the deposited
sequences are 3′-partial (median 79 nt vs the CSV's ~160 nt full Alu), each
AluACA sits *inside* its CSV locus sharing one end; 119 of 122 are fully
contained.

The 222 non-matches are **not** a coordinate-shift artifact — they sit a
median of 1.18 Mb from the nearest CSV locus. Jády's 348 AluACAs and
napRNAdb's 543 Alu/L1 ACA entries are two overlapping but substantially
independent sets, sharing only ~120 loci. Worth knowing before treating
either as the canonical AluACA list.

## The 4 unresolved

- **AluACA96** (APC) — no exact match anywhere in hg38; the deposited
  sequence likely carries a SNP or otherwise differs from the reference.
- **AluACA69**, **AluACA121** — both match `NBPF8` *and* `NBPF9`; the NBPF
  family is highly duplicated and exact matching cannot separate the paralogs.
- **AluACA98** — three identical copies in a chrX tandem array.


## Gotchas encountered (all handled in the scripts)

1. **Getting the accession range right.** `HE` accessions are 2 letters + 6
   digits, and the full AluACA series is **`HE855917`–`HE856264` = 348 records
   = AluACA1–348**, verified at both boundaries (`HE855917` = AluACA1,
   `HE856264` = AluACA348; `HE855915`/`HE855916` are an unrelated
   *Cryptococcus laurentii* rpb1 submission). Starting at `HE855919` silently
   drops AluACA1 and AluACA2 — the count still looks plausible (346), so this
   is easy to miss.
2. **Entrez field qualifiers must be URL-encoded.** Sending `[ACCN]` raw
   returns an empty body with exit code 0 — it looks like a network failure
   but is not. Use `%5BACCN%5D`.
3. **PDF page-break form feeds** (`\f`) prefix the first line of each page, so
   an anchored `/^AluACA/` match silently drops 4 entries.
4. **`AluACA105`** carries a stray trailing `AluACA106,` from the original Word
   layout, and the real 106 entry is typo'd `ALuACA106` (capital L).
5. **66 gene symbols are 2012-era aliases** (`LEPRE1`→`P3H1`, `MLL3`→`KMT2C`,
   `BAT3`→`BAG6`, …), resolved via NCBI Gene and validated against GENCODE.
6. **Duplicate GENCODE entries** for one symbol can yield identical
   coordinates (e.g. `CA5BP1`), which looks like a multi-hit; deduplicated.

## Requirements

- `bedtools` and `python3` (>=3.8, standard library only)
- system `pdftotext` (poppler), for step 02
- an indexed hg38 genome FASTA (`.fa` + `.fai`)
- a GENCODE annotation GTF
- a UCSC RepeatMasker table, for step 09 only
- network access for steps 01 and 04

## Running it elsewhere

Nothing is hard-coded to a particular machine. Every path resolves in three
tiers - built-in default, then environment variable, then command-line flag -
so the same checkout runs on a server without editing any script:

```bash
bash scripts/run_all.sh \
  --proj     /data/aluaca \
  --out      /results/aluaca \
  --work     "$TMPDIR/aluaca" \
  --hg38-fa  /ref/hg38/GRCh38.primary_assembly.genome.fa \
  --gencode  /ref/hg38/gencode.v47.annotation.gtf.gz \
  --rmsk     /ref/hg38/rmsk.tsv.gz
```

Environment variables work equally well, which suits job schedulers:

```bash
export PROJ=/data/aluaca OUT=/results/aluaca WORK=$TMPDIR/aluaca
export HG38_DIR=/ref/hg38 RMSK=/ref/hg38/rmsk.tsv.gz
bash scripts/run_all.sh
```

The flags are parsed by `config.sh`, which every step sources, so they work on
an individual step as well as on `run_all.sh`:

```bash
bash scripts/09_add_repeat_family.sh --out /results/aluaca --rmsk /ref/rmsk.tsv.gz
bash scripts/run_all.sh --help          # full list
```

Details worth knowing:

- **`--proj` vs `--out`.** `--proj` is where the *inputs* live (the PDF, the
  CSV, `deps/`); `--out` is where deliverables are written and defaults to
  `$PROJ`. Split them to keep a read-only input tree. `--work` holds
  intermediates and can point at node-local scratch.
- **Tools.** The `deps/` pixi env is used when present, otherwise `bedtools`
  and `python3` are taken from `PATH` - so a module load or container works
  unchanged. `--bedtools` / `--python` override individually.
- **Reference filenames are guessed.** Given only `--hg38-dir`, the first
  `*.fa`/`*.fasta` and the first `gencode*.gtf.gz` inside it are used, so a
  different GENCODE release needs no edit. `--hg38-fa` / `--gencode` pin them
  explicitly.
- **`--rmsk` accepts gzipped or plain.** `run_all.sh` skips step 09 with a
  notice rather than failing if no RepeatMasker table is configured.
- **Preflight.** `config.sh` checks every required path (including the `.fai`)
  before any step runs and names each missing one with how to set it, so a
  misconfigured run fails in the first second rather than 15 minutes in.

Verified by running the union steps with `--out`, `--work`, `--hg38-dir`,
`--csv` and `--rmsk` all pointed outside the project: byte-identical outputs,
and nothing written into the project directory.

## Union set (`08_build_union.sh`)

Merges the placed Jady AluACAs with the napRNAdb loci, into a single
non-redundant set of 765 intervals.

| file | contents |
|---|---|
| `AluACA_union_nr.bed` | 765 intervals, BED6+1 (with `repeat_family`) |
| `AluACA_union_nr.fasta` | the same 765, stranded sequence from hg38 |
| `AluACA_union_nr.tsv` | BED6 + `source` + `repeat_family` |

**Where the two sets share a locus, the longer interval wins** and both
identifiers are joined in the name. This matters because the deposited AluACA
sequences are 3'-partial (median 79 nt) while the CSV carries the ~160 nt full
element, so the same locus is described at two different extents. The `source`
column records which set each interval's coordinates actually came from:

| `source` | n | coordinates from |
|---|---|---|
| `naprnadb_only` | 422 | CSV, no AluACA at this locus |
| `jady_aluaca` | 222 | AluACA, no CSV locus here |
| `both_naprnadb` | 120 | shared locus, CSV interval was longer |
| `both_jady` | 1 | shared locus, AluACA interval was longer |

Zero residual overlaps on either strand, zero duplicate names.

The 121 shared loci pair 1:1 - no AluACA matches two CSV loci or vice versa.
The CSV interval is longer in 119, AluACA15 wins by 1 bp, and AluACA21 /
hsa-novel-ACA-52 are coordinate-identical (the tie is a no-op).

Two things the numbers hide:

- **344 became 343.** `AluACA88` (HE856004) and `AluACA345` (HE856261) are
  byte-identical sequences deposited twice, both in `RBCK1` at
  chr20:416898-416974 - a duplicate in the original submission, not a mapping
  error. Collapsed to one row.
- **The two sources are not length-comparable.** Jady-only intervals are
  71-162 nt; napRNAdb intervals have a median of 182 nt and run to 6859 nt.
  Anything that scales with interval length - coverage, peak overlap - needs
  normalising or the union will look source-biased.

Verified after the longest-wins step: of the 344 deposited sequences, 224 are
byte-identical to their union record and 118 are exactly contained within a
longer one. The only loss is `AluACA88`/`AluACA345`, whose interval starts
10 bp left of its CSV locus - taking the longer CSV interval (106 nt vs 76 nt)
clips those 10 bp. If that locus matters to you, its full span is
chr20:416898-417014.

### Repeat annotation (`09_add_repeat_family.sh`)

Adds a `repeat_family` column to both union files, making `AluACA_union_nr.bed`
a BED6+1. It holds the RepeatMasker `repName` of every overlapping element,
ordered by overlap length: `AluSx`, or `AluJb|AluSx1` where a locus spans two.

`repName`, not `repFamily` — RepeatMasker calls every Alu's `repFamily` just
"Alu", so the informative J/S/Y subfamily level only exists in `repName`.
Overlap is computed without a strand requirement, and non-Alu elements are
shown rather than hidden.

**All 765 loci sit in an annotated repeat; none are repeat-free.** All 223
Jády-only AluACAs are Alu-dominant, as they must be. The 39 that are not are
napRNAdb-derived and are mostly L1 (the CSV is an Alu *and LINE1* set) plus a
few FLAM/FRAM Alu monomers and MLT/MER LTR fragments. Loci spanning more than
one element are listed in overlap order, e.g. `L1MC4|L1MC4a|AluJb`.

Cross-checked against napRNAdb's own repeat call (CSV column 9) for all 422
shared loci: **zero genuine disagreements.** The CSV's element is the top hit
by overlap in 389, and one of the secondary elements in the other 33 — so the
two annotations differ only in which element of a multi-element locus gets
called dominant, never in identity.
