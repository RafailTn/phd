#!/usr/bin/env bash
# Build the two STAR indices the sno-chimeras pipeline needs.
#
#   1. repeat elements : RepBase human consensus (humrep.ref + humsub.ref)
#   2. genome          : primary assembly + GENCODE splice junctions
#
# Both are used by the --mask_repeat_genome_read stage to discard reads that
# align end-to-end (i.e. non-chimeric reads), and the genome index is reused to
# place the non-guide arm of each candidate chimera.
#
#   bash src/chimeric/build_indices.sh                    # hg38, sparse
#   bash src/chimeric/build_indices.sh --species hg19 --dense
#   bash src/chimeric/build_indices.sh --genome-fa /ref/x.fa --gencode /ref/x.gtf.gz
#
# Build-agnostic: --species selects the output directory and --genome-fa /
# --gencode say what to build it from, so the same script produces the hg38 and
# hg19 indices. Both are idempotent -- an index with an SA file is left alone.
#
# Dense (--dense, sparseD=1) is preferred whenever the machine allows it: the
# sparse array changes where ~6% of short arms are placed, so a sparse index
# cannot be compared against a published dense run without carrying that term in
# the residual. Sparse needs ~17 GB to align, dense ~29 GB.
#
# Run --help for every flag.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

cfg_need "STAR" "$(command -v STAR || true)" "--bin, or put STAR on PATH"
cfg_check

echo "=== [1/2] repeat element index ==="
if [ -f "$REPEAT_INDEX/SA" ]; then
  echo "    exists, skipping"
else
  cfg_need "RepBase FASTA" "$REPEAT_FA" "--repeat-fa"
  cfg_check
  mkdir -p "$REPEAT_INDEX"
  # genomeSAindexNbases = min(14, log2(1.37e6)/2 - 1) ~= 9 for a 1.4 Mbp reference
  STAR --runMode genomeGenerate --runThreadN "$CPUS" \
       --genomeDir "$REPEAT_INDEX" \
       --genomeFastaFiles "$REPEAT_FA" \
       --genomeSAindexNbases 9 \
       --outFileNamePrefix "$REPEAT_INDEX/"
fi

echo "=== [2/2] $SPECIES genome index (sparseD=$SPARSE_D) ==="
if [ -f "$GENOME_INDEX/SA" ]; then
  echo "    exists, skipping"
else
  cfg_need "genome FASTA"    "$GENOME_FA" "--genome-fa, or run fetch_refs.sh $SPECIES"
  cfg_need "GENCODE GTF"     "$GENCODE"   "--gencode, or run fetch_refs.sh $SPECIES"
  cfg_check
  mkdir -p "$GENOME_INDEX"
  # STAR reads a gzipped GTF only via --readFilesCommand, which genomeGenerate
  # does not honour, so decompress both inputs up front if needed.
  FA=$GENOME_FA
  if [[ "$FA" == *.gz ]]; then
    FA=$REF/$(basename "${GENOME_FA%.gz}")
    [ -f "$FA" ] || { echo "    decompressing genome ..."; pigz -dc "$GENOME_FA" > "$FA"; }
  fi
  GTF=$GENCODE
  if [[ "$GTF" == *.gz ]]; then
    GTF=$REF/$(basename "${GENCODE%.gz}")
    [ -f "$GTF" ] || { echo "    decompressing annotation ..."; pigz -dc "$GENCODE" > "$GTF"; }
  fi
  STAR --runMode genomeGenerate --runThreadN "$CPUS" \
       --genomeDir "$GENOME_INDEX" \
       --genomeFastaFiles "$FA" \
       --sjdbGTFfile "$GTF" \
       --sjdbOverhang "$SJDB_OVERHANG" \
       --genomeSAsparseD "$SPARSE_D" \
       --limitGenomeGenerateRAM "$GEN_RAM" \
       --outFileNamePrefix "$GENOME_INDEX/"
fi
echo "=== done ==="
du -sh "$REPEAT_INDEX" "$GENOME_INDEX"
