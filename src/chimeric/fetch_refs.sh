#!/usr/bin/env bash
# Download the references a build needs, into $REF/<build>.
#
#   bash src/chimeric/fetch_refs.sh              # the current --species
#   bash src/chimeric/fetch_refs.sh hg19         # for arm1 / arm3
#   bash src/chimeric/fetch_refs.sh hg19 hg38    # both
#   bash src/chimeric/fetch_refs.sh --ref /scratch/ref hg19
#
# Per build:
#   primary assembly FASTA        ~850-900 MB gz -> ~3.1 GB
#   GENCODE annotation GTF        ~59 MB    (hg19 uses the vNNlift37 mapping, so
#                                            the annotation release is the same
#                                            across builds and only the
#                                            coordinates differ)
#   rmsk.txt.gz -> rmsk.<build>.bed  ~155 MB  RepeatMasker, for the Alu flags
#
# and, for hg19 only, the hg38->hg19 chain used to lift the guide-locus BED.
#
# Which filenames belong to which build is config.sh's business, not this
# script's -- build_indices.sh and run_all.sh read the same variables, so the
# three cannot drift apart.
#
# These are the exact files the workstation run used: the served Content-Length
# for the hg38 GTF (59,075,910) and rmsk.txt.gz (155,633,856) match the local
# copies byte-for-byte, so an hg38 arm run elsewhere differs from the
# workstation's only in the suffix-array density.
#
# Every fetch is skipped if the file is already there, so this is safe to
# re-run and is the second step of run_all.sh.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

cfg_need "curl" "$(command -v curl || true)" "put curl on PATH"
cfg_need "pigz" "$(command -v pigz || true)" "--bin, or put pigz on PATH"
cfg_check

get () {  # get <url> <dest>
  if [ -s "$2" ]; then echo "    have $(basename "$2")"; else
    echo "    fetching $(basename "$2") ..."
    curl -fL --retry 3 -o "$2.part" "$1"
    mv "$2.part" "$2"
  fi
}

fetch_build () {
  local build=$1
  # Re-resolve config.sh's reference table for this build rather than restating
  # it: SPECIES drives GENOME_FA_NAME, GENCODE_NAME and GENCODE_SUBDIR.
  local fa gtf sub names
  names="$(cfg_ref_names "$build")" || {
    echo "unknown build: $build (expected hg19 | hg38)" >&2; return 1; }
  read -r fa gtf sub <<<"$names"
  [ "$sub" = "-" ] && sub=

  local dir=$REF/$build
  mkdir -p "$dir"
  get "$GENCODE_URL/$sub$fa.gz"  "$dir/$fa.gz"
  get "$GENCODE_URL/$sub$gtf"    "$dir/$gtf"
  if [ "$build" = hg19 ]; then
    get "$UCSC_URL/hg38/liftOver/hg38ToHg19.over.chain.gz" "$dir/hg38ToHg19.over.chain.gz"
  fi
  get "$UCSC_URL/$build/database/rmsk.txt.gz" "$dir/rmsk.$build.txt.gz"

  if [ ! -s "$dir/$fa" ]; then
    echo "    decompressing genome ..."
    pigz -dc "$dir/$fa.gz" > "$dir/$fa"
  fi

  # rmsk.txt is the UCSC table dump: bin(1), swScore(2), ..., genoName(6),
  # genoStart(7), genoEnd(8), genoLeft(9), strand(10), repName(11),
  # repClass(12), repFamily(13).
  if [ ! -s "$REF/rmsk.$build.bed" ]; then
    echo "    converting rmsk to BED ..."
    pigz -dc "$dir/rmsk.$build.txt.gz" \
      | awk -v OFS='\t' '{print $6,$7,$8,$11"|"$12"|"$13,0,$10}' \
      | sort -k1,1 -k2,2n > "$REF/rmsk.$build.bed"
  fi

  # The guide-locus BED is in hg38 coordinates; lift it so the false-chimera flag
  # (target landing inside the guide's own locus) still works on hg19.
  if [ "$build" = hg19 ] && [ -s "$REF/guide_loci.hg38.bed" ] && [ ! -s "$REF/guide_loci.hg19.bed" ]; then
    echo "    lifting guide loci hg38 -> hg19 ..."
    liftOver "$REF/guide_loci.hg38.bed" "$dir/hg38ToHg19.over.chain.gz" \
             "$REF/guide_loci.hg19.bed" "$dir/guide_loci.unmapped.bed" || true
    echo "    lifted $(wc -l < "$REF/guide_loci.hg19.bed") of $(wc -l < "$REF/guide_loci.hg38.bed") loci"
  fi
  echo "=== $build references ready in $dir ==="
}

# Bare positionals name the builds; with none, use the configured --species.
BUILDS="$*"
for b in ${BUILDS:-$SPECIES}; do fetch_build "$b"; done
