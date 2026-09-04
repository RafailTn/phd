#!/usr/bin/env bash
# Re-fetch the two runs on the server rather than shipping 1.5 GB over the wire.
#
#   SRR30692552  DKC1 chimeric eCLIP, IP
#   SRR30692553  DKC1 chimeric eCLIP, input control
#
# --seq-defline '@$sn' keeps the original Illumina read names, which the UMI
# handling and every cross-run read-name comparison depend on. Do not drop it.
set -euo pipefail
PROJ=${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
DIR=${DIR:-$PROJ/work/chimeric}
CPUS=${CPUS:-$(nproc)}
[ -d "$PROJ/deps/.pixi/envs/default/bin" ] && export PATH=$PROJ/deps/.pixi/envs/default/bin:$PATH
mkdir -p "$DIR"
for srr in ${SRRS:-SRR30692552 SRR30692553}; do
  if [ -s "$DIR/$srr.fastq.gz" ]; then echo "have $srr"; continue; fi
  echo "=== $srr ==="
  fasterq-dump "$srr" --outdir "$DIR" --temp "$DIR" --threads "$CPUS" \
               --seq-defline '@$sn' --force
  pigz -p "$CPUS" "$DIR/$srr.fastq"
done
ls -la "$DIR"/*.fastq.gz
