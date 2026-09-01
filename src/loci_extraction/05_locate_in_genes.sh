#!/usr/bin/env bash
# Step 05 - Locate each AluACA by exact sequence search INSIDE its annotated
# host gene.
#
# Why not just align to the genome? These are ~76-80 nt fragments of Alu
# consensus; a genome-wide aligner returns hundreds of equally good hits per
# sequence because the genome holds >1M Alu copies. The host-gene assignment
# from Supplemental Table 3 is what makes the mapping unique - so we restrict
# the search space to that gene's span and require an exact match.
#
# Produces work_map/located.tsv with one row per hit and a status column:
#   unique / multi_in_gene / not_found_in_gene / no_gene

source "$(dirname "$0")/config.sh"

echo "[05] building host gene spans ..."
awk -F'\t' 'NR>1 && $2!="NA" && $3!="UNRESOLVED"{print $2}' "$WORK/id_gene.tsv" \
  | sort -u > "$WORK/need_genes.txt"

awk -F'\t' -v OFS='\t' 'NR==FNR{n[$1]=1; next} ($4 in n){print $1,$2,$3,$4,".",$6}' \
  "$WORK/need_genes.txt" "$WORK/genes_hg38.bed" \
  | sort -k1,1 -k2,2n > "$WORK/gene_spans.bed"

echo "  symbols needed: $(wc -l < "$WORK/need_genes.txt")   GENCODE loci: $(wc -l < "$WORK/gene_spans.bed")"
echo "  symbols with no locus: $(comm -23 <(sort "$WORK/need_genes.txt") \
                                          <(cut -f4 "$WORK/gene_spans.bed" | sort -u) | tr '\n' ' ')"
echo "  total span: $(awk '{s+=$3-$2}END{printf "%.1f Mb",s/1e6}' "$WORK/gene_spans.bed")"

echo "[05] extracting gene sequences from hg38 ..."
# No -s: we pull plus-strand sequence and infer strand from which of
# seq / revcomp(seq) matches.
"$BEDTOOLS" getfasta -fi "$HG38_FA" -bed "$WORK/gene_spans.bed" -tab -bedOut \
  > "$WORK/gene_seqs.tsv"
echo "  rows: $(wc -l < "$WORK/gene_seqs.tsv")"

echo "[05] searching ..."
"$PYTHON" "$(dirname "$0")/lib_locate.py"   # inherits WORK/FASTA from config.sh
