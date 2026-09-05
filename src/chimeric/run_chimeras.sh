#!/usr/bin/env bash
# Identify AluACA/snoRNA-guided chimeric reads in DKC1 chimeric eCLIP.
#
# Usage:
#   bash src/chimeric/run_chimeras.sh <SRR> [options] [-- extra sno-chimeras args]
#
#   bash run_chimeras.sh SRR30692552                      # hg38, merged catalogue
#   bash run_chimeras.sh --arm arm1 SRR30692552           # hg19, plain catalogue
#   bash run_chimeras.sh --species hg19 --source plain SRR30692552
#   bash run_chimeras.sh SRR30692552 --outdir /scratch/out --fastq /data/x.fastq.gz
#
# Every path is settable by flag or environment variable; run --help for the
# full list. The legacy positional form `run_chimeras.sh <SRR> [outdir] [fastq]`
# still works.
#
# The guide ("source") catalogue defaults to the merged AluACA+snoRNA union --
# see config.sh for why the real snoRNAs are kept in alongside the AluACAs.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

# Positionals config.sh did not consume: <SRR> [outdir] [fastq]. The last two
# are the legacy form; --outdir / --fastq are the preferred spelling.
SRR=${1:?usage: run_chimeras.sh <SRR> [options] [-- extra sno-chimeras args]}
OUTDIR=${OUTDIR:-${2:-}}
FASTQ=${FASTQ:-${3:-}}
# Legacy 4th-and-beyond positionals still reach sno-chimeras.py, as does
# anything after a literal -- (which config.sh collected into _cfg_extra).
EXTRA=("${@:4}" ${_cfg_extra[@]+"${_cfg_extra[@]}"})

# Results go to $OUT/<arm>/<SRR> when an arm is named, $OUT/<SRR> otherwise.
if [ -z "$OUTDIR" ]; then
  OUTDIR="$OUT${ARM:+/$ARM}/$SRR"
fi
FASTQ="${FASTQ:-$WORK/$SRR.fastq.gz}"

mkdir -p "$OUTDIR"
cfg_need "FASTQ for $SRR"     "$FASTQ"              "--fastq, or run fetch_fastq.sh"
cfg_need "source catalogue"   "$SOURCE_FASTA"       "--source-fasta or --source"
cfg_need "adapter FASTA"      "$ADAPTERS"           "--adapters"
cfg_need "genome STAR index"  "$GENOME_INDEX/SA"    "--genome-index, or run build_indices.sh"
cfg_need "repeat STAR index"  "$REPEAT_INDEX/SA"    "--repeat-index, or run build_indices.sh"
for f in $TARGET_FASTA; do
  cfg_need "target catalogue" "$f"                  "--target-fasta"
done
cfg_check

SOURCE_FASTA=$(readlink -f "$SOURCE_FASTA")
echo "### $SRR | species=$SPECIES${ARM:+ | arm=$ARM} | source=$(basename "$SOURCE_FASTA") ($(grep -c '^>' "$SOURCE_FASTA") records)"

# sno-chimeras.py validates its arguments and only then chdir's into --outdir, so a
# relative FASTQ path passes validation and is dead by the time umi_tools opens it --
# which it does silently, yielding an empty run rather than an error.
FASTQ=$(readlink -f "$FASTQ")
OUTDIR=$(readlink -f "$OUTDIR")

# Opt-in only. sno-chimeras.py uses --source_rna_bed to drop targets landing
# back inside the guide's own locus; the runs behind results/ did not pass it,
# and turning it on silently would change the chimera counts. annotate_chimeras.py
# flags those reads after the fact instead, which is why $GUIDE_BED is not
# wired in here.
SOURCE_BED_ARG=()
if [ -n "${SOURCE_BED:-}" ]; then
  cfg_need "guide-locus BED" "$SOURCE_BED" "--source-bed"
  SOURCE_BED_ARG=(--source_rna_bed "$SOURCE_BED")
fi

"$PYTHON" "$SRC/sno-chimeras.py" "$FASTQ" \
  --uid "$SRR" \
  --strategy chimeric \
  --adapters_fasta "$ADAPTERS" \
  --umi_pattern NNNNNNNNNN \
  --randomer_length 10 \
  --min_length 24 \
  --species "$SPECIES" \
  --genome_star_index "$GENOME_INDEX" \
  --repeat_star_index "$REPEAT_INDEX" \
  --mask_repeat_genome_read \
  --source_rna_fasta "$SOURCE_FASTA" \
  --source_rna_tag snoRNA \
  --target_rna_fasta $TARGET_FASTA \
  --target_rna_tag $TARGET_TAGS \
  --target_rRNA_tag "$TARGET_RRNA_TAG" \
  ${SOURCE_BED_ARG[@]+"${SOURCE_BED_ARG[@]}"} \
  --minimum_chimeras 1000 \
  --cpus "$CPUS" \
  --keep \
  --outdir "$OUTDIR" \
  ${EXTRA[@]+"${EXTRA[@]}"}
