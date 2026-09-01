#!/usr/bin/env bash
# Step 07 - Merge the two localisation passes into a final coordinate set and
# compare it to the napRNAdb CSV.
#
# Confidence tiers kept in the `located_by` column:
#   gene_anchored              exact match inside the annotated host gene
#   genome_unique              single exact match genome-wide
#   genome_disambig_intergenic AluACA11 - two genome hits, and Table 3 says
#                              "No apparent host gene", so the intergenic hit
#                              is taken and the CRY2 hit rejected
#
# Sequences with several equally good genome hits (NBPF paralogs, chrX arrays)
# and the one with no hit at all are written to AluACA_unresolved.tsv instead
# of being forced into the main table.

source "$(dirname "$0")/config.sh"

echo "[07] assembling final coordinates ..."
{
  # gene-anchored; sort -u collapses duplicate GENCODE entries for one symbol
  # that yield identical coordinates (e.g. CA5BP1)
  awk -F'\t' -v OFS='\t' 'NR>1 && $5!="-" && ($11=="unique"||$11=="multi_in_gene"){
        print $5,$6,$7,$1,$3,$8,"gene_anchored"}' "$WORK/located.tsv" | sort -u
  awk -F'\t' -v OFS='\t' 'NR>1 && $2!="-" && $6==1{
        print $2,$3,$4,$1,"-",$5,"genome_unique"}' "$WORK/genome_hits.tsv"
  awk -F'\t' -v OFS='\t' 'NR>1 && $1=="AluACA11" && $2=="chr11" && $3==72150140{
        print $2,$3,$4,$1,"-",$5,"genome_disambig_intergenic"}' "$WORK/genome_hits.tsv"
} | sort -k1,1 -k2,2n > "$WORK/aluaca_hg38.bed"
echo "  located: $(wc -l < "$WORK/aluaca_hg38.bed")"

echo "[07] matching against $(basename "$CSV") ..."
awk -F',' -v OFS='\t' 'NR>1{print $3,$4,$5,$2,$8,$7}' "$CSV" \
  | sort -k1,1 -k2,2n > "$WORK/csv_loci.bed"

"$BEDTOOLS" intersect -a "$WORK/aluaca_hg38.bed" -b "$WORK/csv_loci.bed" -loj \
  | awk -F'\t' -v OFS='\t' '{m=($8=="."?"-":$11); mc=($8=="."?"-":$8":"$9"-"$10);
        print $4,$1,$2,$3,$6,$5,$7,m,mc}' \
  | sort -u > "$WORK/final_body.tsv"

{
  printf "aluaca_id\taccession\ttable3_gene\thg38_chrom\thg38_start\thg38_end\tstrand\tlocated_by\tcsv_napRNA_ID\tcsv_locus\n"
  awk -F'\t' -v OFS='\t' '
    NR==FNR{acc[$2]=$1; next}
    FNR==1 && FILENAME~/id_gene/{next}
    FILENAME~/id_gene/{g[$1]=$2; next}
    {print $1,(acc[$1]?acc[$1]:"-"),(g[$1]?g[$1]:"-"),$2,$3,$4,$5,$7,$8,$9}
  ' <(grep '^>' "$FASTA" | sed -E 's/^>([A-Z0-9.]+) .*(AluACA[0-9]+).*/\1\t\2/') \
    "$WORK/id_gene.tsv" "$WORK/final_body.tsv" \
  | sort -t$'\t' -k1.7n
} > "$OUT/AluACA_hg38_coordinates.tsv"

awk -F'\t' -v OFS='\t' 'NR>1{print $4,$5,$6,$1,"0",$7}' "$OUT/AluACA_hg38_coordinates.tsv" \
  | sort -k1,1 -k2,2n > "$OUT/AluACA_hg38.bed"

# --- sequences that could not be placed unambiguously -------------------
{
  printf "aluaca_id\ttable3_gene\treason\tcandidate_loci\n"
  awk -F'\t' 'NR>1 && $6==0{print $1}' "$WORK/genome_hits.tsv" | while read -r a; do
    g=$(awk -F'\t' -v A="$a" '$1==A{print $2}' "$WORK/id_gene.tsv")
    printf "%s\t%s\tno exact match anywhere in hg38 (likely SNP/assembly diff vs deposited seq)\t-\n" "$a" "$g"
  done
  awk -F'\t' 'NR>1 && $6>1{print $1}' "$WORK/genome_hits.tsv" | sort -u | while read -r a; do
    g=$(awk -F'\t' -v A="$a" '$1==A{print $2}' "$WORK/id_gene.tsv")
    loci=$(awk -F'\t' -v A="$a" 'NR>1 && $1==A{printf "%s:%s-%s(%s);",$2,$3,$4,$5}' "$WORK/genome_hits.tsv")
    [ "$a" = "AluACA11" ] && continue   # disambiguated above
    printf "%s\t%s\tmulti-mapping paralog family; exact match at >1 locus\t%s\n" "$a" "$g" "$loci"
  done
} > "$OUT/AluACA_unresolved.tsv"

echo "[07] summary"
awk -F'\t' 'NR>1{lb[$8]++; if($9!="-")m++}
 END{printf "  rows: %d\n  matched to CSV: %d\n  not in CSV: %d\n",NR-1,m,NR-1-m;
     for(k in lb) printf "  located_by %-28s %s\n",k,lb[k]}' "$OUT/AluACA_hg38_coordinates.tsv"
echo "  CSV loci hit by >=1 AluACA: $("$BEDTOOLS" intersect -a "$WORK/csv_loci.bed" -b "$WORK/aluaca_hg38.bed" -u | wc -l) / $(( $(wc -l < "$CSV") - 1 ))"
echo "[07] wrote $OUT/AluACA_hg38_coordinates.tsv, $OUT/AluACA_hg38.bed, $OUT/AluACA_unresolved.tsv"
