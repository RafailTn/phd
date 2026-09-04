#!/usr/bin/env bash
# Identify AluACA/snoRNA-guided chimeric reads in DKC1 chimeric eCLIP, on hg38.
#
# Usage:  bash src/chimeric/run_chimeras.sh <SRR> [outdir] [fastq]
#
# The source ("guide") catalogue is data/AluACA_snoRNA_merged_nr.fasta -- the 765-locus
# AluACA union merged with the 1936 snoRNA records of snoRNA.txt.fa. Keeping the real
# snoRNAs in alongside the AluACAs is deliberate: bowtie2 scores a read against the whole
# catalogue at once, so an AluACA call means the AluACA beat every snoRNA, rather than
# being the only thing on offer. AluACA records are the ones with a .id3xxx suffix.
set -euo pipefail

PROJ=${PROJ:-/home/rafail/Downloads/phd}
SRR=${1:?usage: run_chimeras.sh <SRR> [outdir] [fastq]}
OUTDIR=${2:-$PROJ/results/chimeric/$SRR}
FASTQ=${3:-$PROJ/work/chimeric/$SRR.fastq.gz}
CPUS=${CPUS:-20}

export PATH=$PROJ/src/chimeric/bin:$PROJ/deps/.pixi/envs/default/bin:$PATH
export CHIMERIC_REF_DIR=$PROJ/ref/chimeric
PY=$PROJ/deps/.pixi/envs/default/bin/python

mkdir -p "$OUTDIR"
[ -f "$FASTQ" ] || { echo "missing FASTQ: $FASTQ" >&2; exit 1; }
# sno-chimeras.py validates its arguments and only then chdir's into --outdir, so a
# relative FASTQ path passes validation and is dead by the time umi_tools opens it --
# which it does silently, yielding an empty run rather than an error.
FASTQ=$(readlink -f "$FASTQ")
OUTDIR=$(readlink -f "$OUTDIR")

"$PY" "$PROJ/src/chimeric/sno-chimeras.py" "$FASTQ" \
  --uid "$SRR" \
  --strategy chimeric \
  --adapters_fasta "$PROJ/ref/chimeric/se.2.round.adapters.fasta" \
  --umi_pattern NNNNNNNNNN \
  --randomer_length 10 \
  --min_length 24 \
  --species hg38 \
  --genome_star_index "$PROJ/ref/chimeric/hg38_star_index" \
  --repeat_star_index "$PROJ/ref/chimeric/repbase_star_index" \
  --mask_repeat_genome_read \
  --source_rna_fasta "$PROJ/data/AluACA_snoRNA_merged_nr.fasta" \
  --source_rna_tag snoRNA \
  --target_rna_fasta "$PROJ/data/rRNA.fa" "$PROJ/data/snRNA.fa" "$PROJ/data/tRNA.fa" \
  --target_rna_tag rRNA snRNA tRNA \
  --target_rRNA_tag rRNA \
  --minimum_chimeras 1000 \
  --cpus "$CPUS" \
  --keep \
  --outdir "$OUTDIR" \
  "${@:4}"
