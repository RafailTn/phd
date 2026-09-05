#!/usr/bin/env bash
# Standalone check (separate from the mapping pipeline):
# do napRNA_Alu_L1_ACA.csv and napRNA_Alu_L1_polyApocketACA.csv overlap?
#
# Result: every one of the 156 polyA-pocket entries has IDENTICAL coordinates
# and strand to an ACA entry, 1:1, and the paired rows agree in every field
# except the ID and the junk "Browser" column. Neither file has any
# within-file overlap. So the ACA file is a strict superset of 543 unique
# loci, and a "keep the longest on overlap" rule never fires - every overlap
# is an exact tie.
#
# Caveat worth remembering: the ACA file has no column marking which 156 of
# its rows are the polyA-pocket ones, so that distinction lives only in the
# pocket file. Recover it any time by joining on Chrom_Start_End.

source "$(dirname "$0")/config.sh"

A="$CSV"
P="$POLYA_CSV"
W="$WORK/csvcheck"; mkdir -p "$W"

echo "=== structure ==="
for f in "$A" "$P"; do
  echo "  $(basename "$f"): rows=$(( $(wc -l < "$f") - 1 ))  fields=$(awk -F',' '{print NF}' "$f" | sort -u | tr '\n' ' ')"
  echo "    Length != End-Start : $(awk -F',' 'NR>1 && ($5-$4)!=$8' "$f" | wc -l)   quoted fields: $(grep -c '"' "$f")"
  echo "    unique IDs=$(awk -F',' 'NR>1{print $2}' "$f" | sort -u | wc -l)  unique coords=$(awk -F',' 'NR>1{print $3"_"$4"_"$5}' "$f" | sort -u | wc -l)"
done

awk -F',' -v OFS='\t' 'NR>1{print $3,$4,$5,$2,$8,$7}' "$A" | sort -k1,1 -k2,2n > "$W/aca.bed"
awk -F',' -v OFS='\t' 'NR>1{print $3,$4,$5,$2,$8,$7}' "$P" | sort -k1,1 -k2,2n > "$W/pocket.bed"

echo "=== between files ==="
echo "  any overlap:         $("$BEDTOOLS" intersect -a "$W/aca.bed" -b "$W/pocket.bed" -wa -wb | wc -l)"
echo "  same strand:         $("$BEDTOOLS" intersect -a "$W/aca.bed" -b "$W/pocket.bed" -wa -wb -s | wc -l)"
echo "  identical coords:    $("$BEDTOOLS" intersect -a "$W/aca.bed" -b "$W/pocket.bed" -wa -wb -f 1.0 -r | wc -l)"
echo "  pocket rows matched: $("$BEDTOOLS" intersect -a "$W/pocket.bed" -b "$W/aca.bed" -u -f 1.0 -r | wc -l)"

echo "=== within files (0 = no self-overlap) ==="
for f in aca pocket; do
  t=$(wc -l < "$W/$f.bed"); n=$("$BEDTOOLS" merge -i "$W/$f.bed" | wc -l)
  echo "  $f: $t intervals -> $n merged ($(( t - n )) overlapping)"
done

echo "=== field-by-field agreement across matched pairs ==="
awk -F',' -v OFS='\t' 'NR>1{print $3"_"$4"_"$5,$0}' "$A" | sort -k1,1 > "$W/a.key"
awk -F',' -v OFS='\t' 'NR>1{print $3"_"$4"_"$5,$0}' "$P" | sort -k1,1 > "$W/p.key"
join -t$'\t' "$W/a.key" "$W/p.key" \
 | awk -F'\t' '{split($2,a,","); split($3,b,","); for(i=1;i<=12;i++) if(a[i]!=b[i]) d[i]++}
   END{split("Browser,napRNAID,Chrom,Start,End,Coverage,Strand,Length,Gene,Region,Conservation,Pubmed",n,",");
       for(i=1;i<=12;i++) printf "  %-12s %s\n",n[i],(d[i]?d[i]" differ":"identical")}'

rm -rf "$W"
