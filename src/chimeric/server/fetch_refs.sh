#!/usr/bin/env bash
# Download the references an arm needs, into $REF/<build>.
#
#   bash fetch_refs.sh hg19          # for arm1 / arm3
#   bash fetch_refs.sh hg38          # for arm0d
#   bash fetch_refs.sh hg19 hg38     # both
#
# Per build:
#   primary assembly FASTA        ~850-900 MB gz -> ~3.1 GB
#   GENCODE v47 annotation GTF    ~59 MB    (hg19 uses the v47lift37 mapping, so
#                                            the annotation release is the same
#                                            across builds and only the
#                                            coordinates differ)
#   rmsk.txt.gz -> rmsk.<build>.bed  ~155 MB  RepeatMasker, for the Alu flags
#
# and, for hg19 only, the hg38->hg19 chain used to lift the guide-locus BED.
#
# These are the exact files the workstation run used: the served Content-Length
# for the hg38 GTF (59,075,910) and rmsk.txt.gz (155,633,856) match the local
# copies byte-for-byte, so an hg38 arm run here differs from the workstation's
# only in the suffix-array density.
set -euo pipefail
PROJ=${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
REF=${REF:-$PROJ/ref/chimeric}
[ -d "$PROJ/deps/.pixi/envs/default/bin" ] && export PATH=$PROJ/deps/.pixi/envs/default/bin:$PATH

G=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47
U=https://hgdownload.soe.ucsc.edu/goldenPath

get () {  # get <url> <dest>
  if [ -s "$2" ]; then echo "    have $(basename "$2")"; else
    echo "    fetching $(basename "$2") ..."
    curl -fL --retry 3 -o "$2.part" "$1"
    mv "$2.part" "$2"
  fi
}

fetch_build () {
  local build=$1 dir=$REF/$1 fa gtf
  mkdir -p "$dir"
  case "$build" in
    hg38) fa=GRCh38.primary_assembly.genome.fa
          get "$G/$fa.gz"                                  "$dir/$fa.gz"
          get "$G/gencode.v47.primary_assembly.annotation.gtf.gz" \
              "$dir/gencode.v47.primary_assembly.annotation.gtf.gz" ;;
    hg19) fa=GRCh37.primary_assembly.genome.fa
          get "$G/GRCh37_mapping/$fa.gz"                   "$dir/$fa.gz"
          get "$G/GRCh37_mapping/gencode.v47lift37.annotation.gtf.gz" \
              "$dir/gencode.v47lift37.annotation.gtf.gz"
          get "$U/hg38/liftOver/hg38ToHg19.over.chain.gz"  "$dir/hg38ToHg19.over.chain.gz" ;;
    *) echo "unknown build: $build (expected hg19 | hg38)" >&2; return 1 ;;
  esac
  get "$U/$build/database/rmsk.txt.gz" "$dir/rmsk.$build.txt.gz"

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

for b in "${@:-hg19}"; do fetch_build "$b"; done
