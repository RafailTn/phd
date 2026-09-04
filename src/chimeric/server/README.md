# Running the dense-index arms on a server

The workstation has 31 GB of RAM. A dense STAR index for GRCh37/38 with splice
junctions needs roughly 32 GB to build and ~29 GB resident to align, so every
arm run locally used `--genomeSAsparseD 2`. That halves the suffix array and
changes where about 6% of short target arms are placed, which is a term we
cannot separate from a genome-build effect. These scripts run the remaining arms
on a machine that can hold a dense index.

**Requirements:** >= 64 GB RAM, and outbound HTTPS to EBI, UCSC and NCBI. Disk:
about 70 GB for the hg19 arms alone, or 150 GB for all three, made up of ~28 GB
per dense index, ~4 GB per genome FASTA (gz plus decompressed), and ~7 GB of kept
intermediates per sample-arm. Drop `--keep` from `run_chimeras.sh` if that last
part is a problem, but then `compare_to_published.py` can no longer attribute a
missing read to a stage.

## The matrix

|          | plain `snoRNA.txt.fa` (1,951) | merged AluACA+snoRNA (2,701)     |
|----------|-------------------------------|----------------------------------|
| **hg38** | arm2 — done locally, sparse   | arm0 — done locally, sparse      |
|          |                               | **arm0d** — same, dense          |
| **hg19** | **arm1** — reproduction       | **arm3** — build stability       |

`arm1` is the cell to judge against the 45,810 genomic chimeras in the published
`DKC1_IP.snoRNA.hg19.chimeras.csv`: same build, same catalogue, dense index.
`arm3` answers whether the AluACA results hold on the authors' build.
`arm0d` is the workstation's own run again, densely: same build, same catalogue,
only the suffix array differs. It is what says whether any headline number in
`RESULTS.md` moves when the sparse index goes away, and it is the arm to run if
you only run one, because every other cell inherits that uncertainty.

The references `fetch_refs.sh hg38` downloads are the same files the workstation
used — the served `Content-Length` for the GENCODE v47 GTF (59,075,910 bytes) and
for UCSC `rmsk.txt.gz` (155,633,856) match the local copies exactly — so arm0d
against arm0 is a clean one-variable comparison.

Two deviations remain in every arm and cannot be closed here: the RepBase
masking index is public RepBase 31.08 rather than the lab's own copy, and the
annotation is GENCODE v47 (lifted to GRCh37 for the hg19 arms) rather than
whatever release the authors used.

## Steps

    # 1. on the workstation
    bash src/chimeric/server/make_bundle.sh
    scp work/chimeric/chimeric-bundle.tar.gz server:

    # 2. on the server
    mkdir phd && tar -xzf chimeric-bundle.tar.gz -C phd && cd phd
    pixi install --manifest-path deps/pixi.toml

    export PROJ=$PWD
    bash src/chimeric/server/fetch_fastq.sh          # ~1.5 GB, both runs
    bash src/chimeric/server/fetch_refs.sh hg19      # genome, GTF, rmsk, liftOver chain
    bash src/chimeric/server/run_arms.sh arm1 arm3

    # and, to re-run the workstation's own configuration with a dense index:
    bash src/chimeric/server/fetch_refs.sh hg38
    bash src/chimeric/server/run_arms.sh arm0d

The bundle deliberately ships neither the FASTQs nor any STAR index: the reads
are re-fetched from SRA and the indices rebuilt, which is the point.

`run_arms.sh` honours `PROJ`, `CPUS`, `SPARSE_D` (leave at 1), `GEN_RAM`, `IP`
and `INPUT`. It builds the index once per genome and skips it if `SA` already
exists, so `arm1 arm3` share one build.

## Timing, extrapolated from the workstation

| step | workstation (20 cores, sparse) | expect on a server |
|---|---|---|
| index build | 17 min | 20–40 min dense |
| IP sample run | 11 min | similar, I/O bound |
| input sample run | 2.5 min | similar |

So `arm1 arm3` is on the order of an hour plus downloads.

## What comes back

Per arm, under `results/chimeric/<arm>/`:

- `<SRR>/<SRR>.snoRNA.<build>.chimeras.csv` — the genomic chimeras
- `<arm>.vs_published.txt` — read-name overlap with the published set, plus a
  stage-by-stage attribution of every published chimera the arm did not
  reproduce (written by `compare_to_published.py`)
- `<SRR>.annotated.tsv` — arm3/arm0d only; guide class and genomic annotation

Copy `results/chimeric/arm*` back to the workstation and the numbers can be
folded into `RESULTS.md`.

## Files

- `fetch_fastq.sh` — re-download both SRA runs, preserving Illumina read names
- `fetch_refs.sh <build>...` — primary assembly, GENCODE v47 (v47lift37 for
  hg19), RepeatMasker, and for hg19 the hg38→hg19 lift of the guide-locus BED
- `run_arms.sh` — build, run, compare, annotate
- `make_bundle.sh` — pack the workstation-side inputs
