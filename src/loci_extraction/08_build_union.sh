#!/usr/bin/env bash
# Step 08 - Non-redundant union of the two AluACA locus sets.
#
#   set A: the 344 placed Jady AluACAs           (AluACA_hg38.bed)
#   set B: the 543 napRNAdb Alu/L1 ACA loci      (napRNA_Alu_L1_ACA.csv)
#
# 121 loci are shared. Because the deposited AluACA sequences are 3'-partial
# (median 79 nt) while the CSV carries the ~160 nt full element, the two sets
# describe the same locus at different extents. Where they overlap the union
# keeps the LONGER interval and joins both identifiers, so no locus is
# represented by a truncated interval.
#
# One genuine within-set duplicate is collapsed first: AluACA88 (HE856004) and
# AluACA345 (HE856261) are byte-identical sequences deposited twice, both in
# RBCK1 at chr20:416898-416974. Kept once, both IDs in the name.
#
# The `source` column records which set each interval's coordinates came from:
#   jady_aluaca     AluACA only, coordinates from AluACA_hg38.bed
#   naprnadb_only   CSV only, coordinates from the CSV
#   both_jady       shared locus, the AluACA interval was the longer one
#   both_naprnadb   shared locus, the CSV interval was the longer one
#
# Writes AluACA_union_nr.bed / .fasta / .tsv to the project root.

source "$(dirname "$0")/config.sh"
W="$WORK/union"; mkdir -p "$W"

awk -F',' -v OFS='\t' 'NR>1{print $3,$4,$5,$2,$8,$7}' "$CSV" | sort -k1,1 -k2,2n > "$W/csv.bed"
sort -k1,1 -k2,2n "$OUT/AluACA_hg38.bed" > "$W/aluaca_raw.bed"

# collapse identical-interval AluACA duplicates, joining the IDs with "|"
"$BEDTOOLS" merge -i "$W/aluaca_raw.bed" -s -c 4,6 -o distinct,distinct \
  | awk -v OFS='\t' '{gsub(",","|",$4); print $1,$2,$3,$4,"0",$5}' \
  | sort -k1,1 -k2,2n > "$W/aluaca.bed"

# shared loci: pair them up, then keep whichever interval is longer
"$BEDTOOLS" intersect -a "$W/aluaca.bed" -b "$W/csv.bed" -wa -wb -s > "$W/pairs.tsv"
awk -F'\t' -v OFS='\t' '
  { la=$3-$2; lc=$9-$8
    if (la > lc) print $1,$2,$3,$4"|"$10,0,$6,"both_jady"
    else         print $7,$8,$9,$4"|"$10,0,$12,"both_naprnadb" }' "$W/pairs.tsv" > "$W/shared.bed7"

# the rest of each set, untouched
"$BEDTOOLS" intersect -a "$W/aluaca.bed" -b "$W/csv.bed" -v -s \
  | awk -v OFS='\t' '{print $1,$2,$3,$4,0,$6,"jady_aluaca"}'   > "$W/jady_only.bed7"
"$BEDTOOLS" intersect -a "$W/csv.bed" -b "$W/aluaca.bed" -v -s \
  | awk -v OFS='\t' '{print $1,$2,$3,$4,0,$6,"naprnadb_only"}' > "$W/csv_only.bed7"

cat "$W/shared.bed7" "$W/jady_only.bed7" "$W/csv_only.bed7" \
  | sort -k1,1 -k2,2n > "$W/union.bed7"

cut -f1-6 "$W/union.bed7" > "$OUT/AluACA_union_nr.bed"
{ printf "chrom\tstart\tend\tname\tscore\tstrand\tsource\n"; cat "$W/union.bed7"; } \
  > "$OUT/AluACA_union_nr.tsv"

# stranded sequence straight from hg38, so both sources are on equal footing
"$BEDTOOLS" getfasta -fi "$HG38_FA" -bed "$OUT/AluACA_union_nr.bed" -s -name \
  | sed -E 's/^>([^:]+)::/>\1 /' > "$OUT/AluACA_union_nr.fasta"

echo "[08] union summary"
cut -f7 "$W/union.bed7" | sort | uniq -c | awk '{printf "  %-16s %s\n",$2,$1}'
echo "  union intervals:               $(wc -l < "$OUT/AluACA_union_nr.bed")"
echo "  fasta records:                 $(grep -c '^>' "$OUT/AluACA_union_nr.fasta")"
echo "  residual same-strand overlaps: $(( $(wc -l < "$OUT/AluACA_union_nr.bed") - $("$BEDTOOLS" merge -s -i "$OUT/AluACA_union_nr.bed" | wc -l) ))"
echo "  duplicate names:               $(( $(cut -f4 "$OUT/AluACA_union_nr.bed" | wc -l) - $(cut -f4 "$OUT/AluACA_union_nr.bed" | sort -u | wc -l) ))"
echo "  AluACA sequence not covered by the chosen interval:"
awk -F'\t' '($9-$8)>=($3-$2) && ($2<$8 || $3>$9) {printf "    %-22s aluaca %s:%d-%d  kept %s:%d-%d\n",$4,$1,$2,$3,$7,$8,$9}' "$W/pairs.tsv"
