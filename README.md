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
