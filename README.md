# PhD analyses

Three pipelines, on AluACA RNAs — a class of Alu-derived H/ACA-like small RNAs — and
whether they behave like real DKC1 guides.

| directory | what it does | run it with |
|---|---|---|
| [`src/loci_extraction/`](src/loci_extraction/) | places the published AluACA sequences on hg38 and builds the 765-locus union catalogue | `bash src/loci_extraction/run_all.sh` |
| [`src/analysis/`](src/analysis/) | compares that catalogue against snoDB's snoRNAs and merges the two into one non-redundant guide set | `python3 src/analysis/collapse_duplicates.py` |
| [`src/chimeric/`](src/chimeric/) | finds AluACA-guided chimeric reads in DKC1 chimeric eCLIP and says what the other arm is | `bash src/chimeric/run_all.sh` |

They run in that order: `loci_extraction` produces the union catalogue, `analysis` merges
it with the snoRNAs, and `chimeric` uses the merged catalogue as its guide set. Each
directory has its own README with the method and the reasoning.

The result is [`results/chimeric/RESULTS.md`](results/chimeric/RESULTS.md), regenerated
from the pipeline output by `make_report.py` rather than edited by hand.

## Conventions

Every path in every pipeline is settable three ways, in increasing precedence: a built-in
default, an environment variable, then a command-line flag. Run any script with `--help`
for its full list. Nothing is hardcoded to one machine — the project root is found by
walking up from the script — so a checkout runs wherever it is put:

```bash
bash src/chimeric/run_all.sh --proj /data/phd --out /results --cpus 32
bash src/loci_extraction/run_all.sh --hg38-dir /ref/hg38 --work "$TMPDIR/aluaca"
```

`src/loci_extraction/config.sh` and `src/chimeric/config.sh` hold each pipeline's
contract; `src/paths.py` is the Python equivalent, shared by the scripts in
`src/analysis/` and `src/chimeric/`.

### Where inputs are found

You do not normally have to say. An input is looked for by name in the conventional
places first, then **anywhere under the project root**, so it is enough to drop a file
somewhere in the repo:

| pipeline | looked for in, in order |
| --- | --- |
| `chimeric` | `data/`, then the project root, then anywhere under it |
| `loci_extraction` | the project root, then `data/`, then anywhere under it |
| `analysis` | the project root, then `data/`, then anywhere under it |

The search **never leaves the project**. The current directory is deliberately not
consulted, so running from elsewhere cannot silently pick up a same-named file; `deps/`,
`.git/`, `__pycache__/` and the STAR indices are skipped. Of several matches the
shallowest wins, ties broken alphabetically, so the choice never depends on filesystem
order. To use a file from outside the repo, name it with its flag.

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

A missing input fails in the first second with a message naming the flag that fixes it,
rather than part-way through a run.

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

The scripts prefer `deps/.pixi/envs/default/bin` when it exists and fall back to `PATH`
otherwise, so a module system, conda environment or container works too.
