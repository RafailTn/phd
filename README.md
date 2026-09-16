# PhD analyses

Three pipelines. They study AluACA RNAs, a class of Alu-derived H/ACA-like small RNAs,
and test whether these RNAs act as real DKC1 guides.

| directory | what it does | run it with |
|---|---|---|
| [`src/loci_extraction/`](src/loci_extraction/) | places the published AluACA sequences on hg38 and builds the 765-locus union catalogue | `bash src/loci_extraction/run_all.sh` |
| [`src/analysis/`](src/analysis/) | compares that catalogue against snoDB's snoRNAs and merges the two into one non-redundant guide set | `python3 src/analysis/collapse_duplicates.py` |
| [`src/chimeric/`](src/chimeric/) | finds AluACA-guided chimeric reads in DKC1 chimeric eCLIP and says what the other arm is | `bash src/chimeric/run_all.sh` |

Run them in that order:

1. `loci_extraction` produces the union catalogue.
2. `analysis` merges the catalogue with the snoRNAs.
3. `chimeric` uses the merged catalogue as its guide set.

Each directory has its own README with the method and the reasoning.

The result is [`results/chimeric_chr25/RESULTS.hg19.md`](results/chimeric_chr25/RESULTS.hg19.md),
on hg19 with the 25-contig index. `make_report.py` writes it from the pipeline's output.
Nobody edits it by hand. The hg38 reports are kept for the record but superseded: GRCh38's
rDNA models make that build unable to reproduce the published calls.

## Conventions

Every path in every pipeline is settable three ways, in increasing precedence: a
built-in default, an environment variable, then a command-line flag. Run any script with
`--help` for its full list. No path is specific to one machine. Each script finds the
project root by walking up from its own location. A checkout therefore runs in any
directory:

```bash
bash src/chimeric/run_all.sh --proj /data/phd --out /results --cpus 32
bash src/loci_extraction/run_all.sh --hg38-dir /ref/hg38 --work "$TMPDIR/aluaca"
```

`src/loci_extraction/config.sh` and `src/chimeric/config.sh` hold each pipeline's
contract. `src/paths.py` is the Python equivalent. The scripts in `src/analysis/` and
`src/chimeric/` share it.

### How the scripts find inputs

You do not normally have to say. Each script looks for an input by name in the
conventional places first, then **anywhere under the project root**. Put the file
anywhere in the repo and the script finds it:

| pipeline | search order |
| --- | --- |
| `chimeric` | `data/`, then the project root, then anywhere under it |
| `loci_extraction` | the project root, then `data/`, then anywhere under it |
| `analysis` | the project root, then `data/`, then anywhere under it |

The search **never leaves the project**. The search ignores the current directory on
purpose, so a run from another directory cannot pick up a same-named file. The search
also ignores `deps/`, `.git/`, `__pycache__/` and the STAR indices. If several files
match, the shallowest one wins. The script breaks a tie alphabetically. The choice never
depends on filesystem order. To use a file from outside the repo, name it with its flag.

### Overriding an input

Each input has a flag, and an environment variable of the same name in upper snake case:

**`src/chimeric/`** — `bash src/chimeric/run_all.sh --help`

| flag | what |
| --- | --- |
| `--source-fasta` / `--source` | the guide catalogue: a path, or `plain` / `merged` |
| `--fastq` | the reads for one run (otherwise `work/chimeric/<SRR>.fastq.gz`) |
| `--alu-fasta` | FASTA whose headers name the AluACA records |
| `--target-fasta` / `--target-tag` | target RNA catalogues, repeatable and paired |
| `--adapters` | second-round adapter FASTA |
| `--published` | published hg19 chimeras CSV, for the comparison |
| `--genome-fa` `--gencode` `--rmsk` `--guide-bed` `--repeat-fa` | references |

**`src/loci_extraction/`** — `bash src/loci_extraction/run_all.sh --help`

| flag | what |
| --- | --- |
| `--csv` | napRNAdb Alu/L1 ACA CSV |
| `--polya-csv` | napRNAdb Alu/L1 polyA-pocket ACA CSV |
| `--pdf` | Jady et al. supplemental PDF |
| `--fasta` | deposited-sequence FASTA (written by step 01) |
| `--snodb` | snoDB catalogue TSV |
| `--hg38-fa` `--gencode` `--rmsk` | references |

**`src/analysis/`** — `python3 src/analysis/<script>.py --help`

| script | input flags |
| --- | --- |
| `collapse_duplicates.py` | `--sno` `--union` `--bed` |
| `snorna_locate.py` | `--fasta` `--genome` `--bed` |
| `snorna_overlap.py` | `--a` `--b` |

```bash
# all equivalent ways to point at one input
bash src/chimeric/run_chimeras.sh --source-fasta /ref/custom.fa SRR30692552
SOURCE_FASTA=/ref/custom.fa bash src/chimeric/run_chimeras.sh SRR30692552
cp /ref/custom.fa data/ && bash src/chimeric/run_chimeras.sh --source-fasta data/custom.fa SRR30692552
```

A script stops in the first second if an input is missing. The message names the flag
that supplies the input. The script does not stop part-way through a run.

## Layout

```
data/     inputs -- catalogues, published tables, SRA runs   (gitignored)
ref/      references -- genomes, annotations, STAR indices   (gitignored)
work/     scratch, FASTQs, logs                              (gitignored)
results/  outputs; only RESULTS.md and the arm comparisons are committed
deps/     the pixi environment that provides every tool
src/      the pipelines
```

## Dependencies

Everything comes from one [pixi](https://pixi.sh) environment:

```bash
pixi install --manifest-path deps/pixi.toml
```

The scripts use `deps/.pixi/envs/default/bin` when that directory exists. If it does not
exist, they use `PATH`. A module system, conda environment or container therefore also
works.
