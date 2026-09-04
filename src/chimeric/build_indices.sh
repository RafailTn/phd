#!/usr/bin/env bash
# Build the two STAR indices the sno-chimeras pipeline needs, for hg38.
#
#   1. repeat elements : RepBase human consensus (humrep.ref + humsub.ref)
#   2. genome          : GRCh38 primary assembly + GENCODE v47 splice junctions
#
# Both are used by the --mask_repeat_genome_read stage to discard reads that
# align end-to-end (i.e. non-chimeric reads), and the genome index is reused to
# place the non-guide arm of each candidate chimera.
set -euo pipefail

PROJ=${PROJ:-/home/rafail/Downloads/phd}
REF=$PROJ/ref/chimeric
HG38_DIR=${HG38_DIR:-/home/rafail/Downloads/hg38}
HG38_FA=${HG38_FA:-$HG38_DIR/GRCh38.primary_assembly.genome.fa}
GENCODE=${GENCODE:-$HG38_DIR/gencode.v47.primary_assembly.annotation.gtf.gz}
CPUS=${CPUS:-20}
# sjdbOverhang = max read length - 1. Reads are 150 nt minus a 10 nt UMI, so 139.
SJDB_OVERHANG=${SJDB_OVERHANG:-139}

export PATH=$PROJ/deps/.pixi/envs/default/bin:$PATH

echo "=== [1/2] repeat element index ==="
if [ -f "$REF/repbase_star_index/SA" ]; then
  echo "    exists, skipping"
else
  mkdir -p "$REF/repbase_star_index"
  # genomeSAindexNbases = min(14, log2(1.37e6)/2 - 1) ~= 9 for a 1.4 Mbp reference
  STAR --runMode genomeGenerate --runThreadN "$CPUS" \
       --genomeDir "$REF/repbase_star_index" \
       --genomeFastaFiles "$REF/repbase/human_repbase.fa" \
       --genomeSAindexNbases 9 \
       --outFileNamePrefix "$REF/repbase_star_index/"
fi

echo "=== [2/2] hg38 genome index (this takes ~1h) ==="
if [ -f "$REF/hg38_star_index/SA" ]; then
  echo "    exists, skipping"
else
  mkdir -p "$REF/hg38_star_index"
  GTF=$GENCODE
  if [[ "$GENCODE" == *.gz ]]; then
    GTF=$REF/$(basename "${GENCODE%.gz}")
    [ -f "$GTF" ] || { echo "    decompressing annotation ..."; pigz -dc "$GENCODE" > "$GTF"; }
  fi
  # genomeSAsparseD 2 halves the suffix-array memory (build and align); needed to
  # fit GRCh38 + sjdb into the 31 GB on this machine.
  STAR --runMode genomeGenerate --runThreadN "$CPUS" \
       --genomeDir "$REF/hg38_star_index" \
       --genomeFastaFiles "$HG38_FA" \
       --sjdbGTFfile "$GTF" \
       --sjdbOverhang "$SJDB_OVERHANG" \
       --genomeSAsparseD 2 \
       --limitGenomeGenerateRAM 26000000000 \
       --outFileNamePrefix "$REF/hg38_star_index/"
fi
echo "=== done ==="
du -sh "$REF"/*_star_index
