#!/usr/bin/env bash
# Download the sequencing runs into $WORK, where run_chimeras.sh looks for them.
#
#   bash src/chimeric/fetch_fastq.sh                       # the IP and its control
#   bash src/chimeric/fetch_fastq.sh --srr SRR30692552     # just one
#   bash src/chimeric/fetch_fastq.sh --work /scratch/fastq
#
#   SRR30692552  DKC1 chimeric eCLIP, IP
#   SRR30692553  DKC1 chimeric eCLIP, input control
#
# Existing files are left alone, so this is safe to re-run and is the first step
# of run_all.sh.
#
# --seq-defline '@$sn' keeps the original Illumina read names, which the UMI
# handling and every cross-run read-name comparison depend on. Do not drop it.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

cfg_need "fasterq-dump" "$(command -v fasterq-dump || true)" "--bin, or put sra-tools on PATH"
cfg_need "pigz"         "$(command -v pigz || true)"         "--bin, or put pigz on PATH"
cfg_check

# Accessions come from bare positionals, from --srr, or from the default pair.
ACCESSIONS="$*"
ACCESSIONS="${ACCESSIONS:-$SRRS}"

mkdir -p "$WORK"
for srr in $ACCESSIONS; do
  if [ -s "$WORK/$srr.fastq.gz" ]; then echo "have $srr"; continue; fi
  echo "=== $srr ==="
  fasterq-dump "$srr" --outdir "$WORK" --temp "$WORK" --threads "$CPUS" \
               --seq-defline '@$sn' --force
  pigz -p "$CPUS" "$WORK/$srr.fastq"
done
ls -la "$WORK"/*.fastq.gz
