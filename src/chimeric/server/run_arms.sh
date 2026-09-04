#!/usr/bin/env bash
# Run the reproduction-matrix arms that need a dense STAR index, on a machine
# with enough RAM for one (>= 64 GB recommended; STAR holds ~29 GB resident for a
# dense GRCh37/38 index with sjdb, and genomeGenerate peaks higher).
#
#   bash run_arms.sh arm1            # hg19 + plain 1951-record snoRNA.txt.fa
#   bash run_arms.sh arm3            # hg19 + merged 2701-record catalogue
#   bash run_arms.sh arm1 arm3       # both, reusing the same index
#   bash run_arms.sh arm0d           # the workstation run again, dense
#
# The 2x2 this completes:
#
#            | plain snoRNA.txt.fa      | merged AluACA+snoRNA
#   ---------+--------------------------+--------------------------
#   hg38     | arm2 (done, sparse)      | arm0 (done, sparse)
#   hg19     | arm1  <- reproduction    | arm3  <- build-stability
#
# arm1 is the number to judge against the published 45,810 genomic chimeras in
# DKC1_IP.snoRNA.hg19.chimeras.csv. arm3 says whether the AluACA findings survive
# a change of build. arm0d re-runs the workstation's own configuration densely and
# so measures the sparse-vs-dense term that every workstation arm carries; the
# references it downloads are byte-identical to the ones used there, so the index
# density is the only thing that differs.
set -euo pipefail

PROJ=${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
REF=${REF:-$PROJ/ref/chimeric}
CPUS=${CPUS:-$(nproc)}
# 1 = dense. This is the whole point of running here rather than on the workstation.
SPARSE_D=${SPARSE_D:-1}
GEN_RAM=${GEN_RAM:-64000000000}
IP=${IP:-SRR30692552}
INPUT=${INPUT:-SRR30692553}
export PROJ REF CPUS

[ -d "$PROJ/deps/.pixi/envs/default/bin" ] && export PATH=$PROJ/deps/.pixi/envs/default/bin:$PATH

PLAIN=$PROJ/data/snoRNA.txt.fa
MERGED=$PROJ/data/AluACA_snoRNA_merged_nr.fasta

build () {  # build <species>
  SPECIES=$1 SPARSE_D=$SPARSE_D GEN_RAM=$GEN_RAM CPUS=$CPUS \
    bash "$PROJ/src/chimeric/build_indices.sh"
}

one () {  # one <arm> <species> <source fasta> <sample> ...
  local arm=$1 species=$2 src=$3; shift 3
  for srr in "$@"; do
    local out=$PROJ/results/chimeric/${arm}/${srr}
    echo "######## $arm $srr start $(date) ########"
    SPECIES=$species SOURCE_FASTA=$src CPUS=$CPUS \
      bash "$PROJ/src/chimeric/run_chimeras.sh" "$srr" "$out"
    echo "######## $arm $srr done $(date) ########"
  done
}

compare () {  # compare <arm> <species> <label>
  local arm=$1 species=$2 label=$3
  echo "######## compare $arm vs published $(date) ########"
  python "$PROJ/src/chimeric/compare_to_published.py" \
    --outdir "$PROJ/results/chimeric/${arm}/${IP}" --uid "$IP" --gtag "$species" \
    --published "$PROJ/data/DKC1_IP.snoRNA.hg19.chimeras.csv" \
    --label "$label" \
    --out "$PROJ/results/chimeric/${arm}.vs_published.txt"
}

annotate () {  # annotate <arm> <species> <sample> ...
  local arm=$1 species=$2; shift 2
  local gtf=$REF/hg19/gencode.v47lift37.annotation.gtf.gz
  local bed=$REF/guide_loci.hg19.bed
  local rmsk=$REF/rmsk.hg19.bed
  if [ "$species" = hg38 ]; then
    gtf=$REF/hg38/gencode.v47.primary_assembly.annotation.gtf.gz
    [ -f "$gtf" ] || gtf=/home/rafail/Downloads/hg38/gencode.v47.primary_assembly.annotation.gtf.gz
    bed=$REF/guide_loci.hg38.bed; rmsk=$REF/rmsk.hg38.bed
  fi
  for srr in "$@"; do
    echo "######## annotate $arm $srr $(date) ########"
    python "$PROJ/src/chimeric/annotate_chimeras.py" \
      --outdir "$PROJ/results/chimeric/${arm}/${srr}" \
      --uid "$srr" --gtag "$species" --tags "rRNA,snRNA,tRNA,$species" \
      --gtf "$gtf" --alu-fasta "$PROJ/data/AluACA_union_nr.fasta" \
      --source-bed "$bed" --rmsk "$rmsk" \
      --out "$PROJ/results/chimeric/${arm}/${srr}.annotated.tsv"
  done
}

for arm in "$@"; do
  case "$arm" in
    arm1)  build hg19; one arm1_hg19_plain  hg19 "$PLAIN"  "$IP"
           compare arm1_hg19_plain  hg19 "arm1: hg19 + plain (dense)" ;;
    arm3)  build hg19; one arm3_hg19_merged hg19 "$MERGED" "$IP" "$INPUT"
           compare arm3_hg19_merged hg19 "arm3: hg19 + merged (dense)"
           annotate arm3_hg19_merged hg19 "$IP" "$INPUT" ;;
    arm0d) build hg38; one arm0d_hg38_merged_dense hg38 "$MERGED" "$IP" "$INPUT"
           compare arm0d_hg38_merged_dense hg38 "arm0d: hg38 + merged (dense)"
           annotate arm0d_hg38_merged_dense hg38 "$IP" "$INPUT" ;;
    *) echo "unknown arm: $arm (expected arm1 | arm3 | arm0d)" >&2; exit 1 ;;
  esac
done
echo ARMS_COMPLETE
