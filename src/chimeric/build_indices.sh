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
# Build-agnostic: SPECIES selects the output directory and GENOME_FA/GENCODE say
# what to build it from, so the same script produces the hg38 and hg19 indices.
#
#   SPECIES     hg38 | hg19        (default hg38)
#   GENOME_FA   primary assembly FASTA (plain or .gz)
#   GENCODE     annotation GTF (plain or .gz)
#   SPARSE_D    STAR --genomeSAsparseD. 1 = dense (needs ~32 GB free RAM to
#               build and ~29 GB to align); 2 = halved suffix array, ~17 GB.
#               Default 2, which is what the 31 GB workstation requires.
#   GEN_RAM     --limitGenomeGenerateRAM in bytes
#   CPUS        threads
#
# Dense is preferred whenever the machine allows it: the sparse array changes
# where ~6% of short arms are placed, so a sparse index cannot be compared
# against a published dense run without carrying that term in the residual.
set -euo pipefail

PROJ=${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}
REF=${REF:-$PROJ/ref/chimeric}
SPECIES=${SPECIES:-hg38}
CPUS=${CPUS:-20}
SPARSE_D=${SPARSE_D:-2}
GEN_RAM=${GEN_RAM:-26000000000}
# sjdbOverhang = max read length - 1. Reads are 150 nt minus a 10 nt UMI, so 139.
SJDB_OVERHANG=${SJDB_OVERHANG:-139}

# Defaults resolve to $REF/<species>, where server/fetch_refs.sh puts things. The
# workstation keeps its hg38 copy outside the repo, so that location wins if it is
# present -- the two are the same files (GENCODE release 47, primary assembly).
WS_HG38=${WS_HG38:-/home/rafail/Downloads/hg38}
case "$SPECIES" in
  hg38) if [ -z "${GENOME_FA:-}" ] && [ -f "$WS_HG38/GRCh38.primary_assembly.genome.fa" ]; then
          GENOME_FA=$WS_HG38/GRCh38.primary_assembly.genome.fa
          GENCODE=${GENCODE:-$WS_HG38/gencode.v47.primary_assembly.annotation.gtf.gz}
        fi
        GENOME_FA=${GENOME_FA:-$REF/hg38/GRCh38.primary_assembly.genome.fa}
        GENCODE=${GENCODE:-$REF/hg38/gencode.v47.primary_assembly.annotation.gtf.gz} ;;
  hg19) GENOME_FA=${GENOME_FA:-$REF/hg19/GRCh37.primary_assembly.genome.fa}
        GENCODE=${GENCODE:-$REF/hg19/gencode.v47lift37.annotation.gtf.gz} ;;
  *)    : "${GENOME_FA:?set GENOME_FA for SPECIES=$SPECIES}"
        : "${GENCODE:?set GENCODE for SPECIES=$SPECIES}" ;;
esac

GENOME_INDEX=${GENOME_INDEX:-$REF/${SPECIES}_star_index}
REPEAT_INDEX=${REPEAT_INDEX:-$REF/repbase_star_index}
REPEAT_FA=${REPEAT_FA:-$REF/repbase/human_repbase.fa}

[ -d "$PROJ/deps/.pixi/envs/default/bin" ] && export PATH=$PROJ/deps/.pixi/envs/default/bin:$PATH
command -v STAR >/dev/null || { echo "STAR not on PATH" >&2; exit 1; }

echo "=== [1/2] repeat element index ==="
if [ -f "$REPEAT_INDEX/SA" ]; then
  echo "    exists, skipping"
else
  [ -f "$REPEAT_FA" ] || { echo "missing $REPEAT_FA" >&2; exit 1; }
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
  [ -f "$GENOME_FA" ] || { echo "missing $GENOME_FA" >&2; exit 1; }
  [ -f "$GENCODE" ]   || { echo "missing $GENCODE" >&2; exit 1; }
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
