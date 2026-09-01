#!/usr/bin/env bash
# Step 03 - Flatten GENCODE v47 gene records to BED6 + a symbol list.
# Used by step 04 (validating 2012-era symbols) and step 05 (gene spans).

source "$(dirname "$0")/config.sh"

echo "[03] flattening GENCODE gene records ..."
zcat "$GENCODE" \
  | awk -F'\t' -v OFS='\t' '$3=="gene"{
      match($9,/gene_name "[^"]+"/); gn=substr($9,RSTART+11,RLENGTH-12);
      match($9,/gene_type "[^"]+"/); gt=substr($9,RSTART+11,RLENGTH-12);
      print $1,$4-1,$5,gn,gt,$7 }' \
  | grep -P "^chr" \
  | sort -k1,1 -k2,2n > "$WORK/genes_hg38.bed"

cut -f4 "$WORK/genes_hg38.bed" | sort -u > "$WORK/gene_names.txt"
echo "[03] genes: $(wc -l < "$WORK/genes_hg38.bed")  unique symbols: $(wc -l < "$WORK/gene_names.txt")"
