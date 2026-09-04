#!/usr/bin/env bash
# Pack everything the server needs, except the FASTQs and the STAR indices
# (which are rebuilt there -- that is the point of the exercise).
#
#   bash src/chimeric/server/make_bundle.sh [outfile.tar.gz]
#
# On the server:
#   mkdir phd && tar -xzf chimeric-bundle.tar.gz -C phd && cd phd
#   pixi install --manifest-path deps/pixi.toml
#   bash src/chimeric/server/fetch_fastq.sh
#   bash src/chimeric/server/fetch_refs.sh hg19
#   bash src/chimeric/server/run_arms.sh arm1 arm3
set -euo pipefail
PROJ=${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
OUT=${1:-$PROJ/work/chimeric/chimeric-bundle.tar.gz}
cd "$PROJ"

FILES=(
  src/chimeric
  deps/pixi.toml deps/pixi.lock
  data/snoRNA.txt.fa
  data/AluACA_snoRNA_merged_nr.fasta
  data/AluACA_union_nr.fasta
  data/rRNA.fa data/snRNA.fa data/tRNA.fa
  data/DKC1_IP.snoRNA.hg19.chimeras.csv
  ref/chimeric/se.2.round.adapters.fasta
  ref/chimeric/repbase/human_repbase.fa
  ref/chimeric/guide_loci.hg38.bed
)
for f in "${FILES[@]}"; do
  [ -e "$f" ] || { echo "missing: $f" >&2; exit 1; }
done

mkdir -p "$(dirname "$OUT")"
tar -czf "$OUT" --exclude='__pycache__' "${FILES[@]}"
echo "wrote $OUT ($(du -h "$OUT" | cut -f1))"
tar -tzf "$OUT" | wc -l | xargs echo "entries:"
