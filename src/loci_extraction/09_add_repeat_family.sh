#!/usr/bin/env bash
# Step 09 - Annotate the union set with the repeat element(s) each locus sits in.
#
# Adds one column, `repeat_family`, holding the RepeatMasker repName of every
# overlapping element ordered by overlap length, e.g. "AluSx" or "AluJb|AluSx1".
# repName is used rather than repFamily because RepeatMasker's repFamily for
# every Alu is just "Alu" - the informative level (AluJ/AluS/AluY subfamilies)
# lives in repName. Non-Alu overlaps (L1, MIR, ...) are shown too rather than
# hidden, and "none" means no annotated repeat at all.
#
# Overlap is computed without a strand requirement; strand concordance between
# the locus and the element is reported in the summary instead.
#
# Rewrites AluACA_union_nr.bed (BED6+1) and AluACA_union_nr.tsv in place.

source "$(dirname "$0")/config.sh"
# RMSK comes from config.sh (--rmsk / $RMSK)
W="$WORK/union"; mkdir -p "$W"
[ -n "${RMSK:-}" ] && [ -f "$RMSK" ] || {
  echo "[09] missing RepeatMasker table: ${RMSK:-<unset>}   (set it with --rmsk)" >&2; exit 1; }

# gzipped or plain
_cat() { case "$RMSK" in *.gz) zcat "$1";; *) cat "$1";; esac; }

# UCSC rmsk table (tab-separated despite the .gtf name):
#   6 genoName  7 genoStart  8 genoEnd  10 strand  11 repName  12 repClass  13 repFamily
echo "[09] converting RepeatMasker table ..."
_cat "$RMSK" \
  | awk -F'\t' -v OFS='\t' 'NR>1{print $6,$7,$8,$11,$12";"$13,$10}' \
  | sort -k1,1 -k2,2n > "$W/rmsk.bed"
echo "  elements: $(wc -l < "$W/rmsk.bed")"

echo "[09] intersecting ..."
"$BEDTOOLS" intersect -a "$W/union.bed7" -b "$W/rmsk.bed" -wao > "$W/annot_raw.tsv"

# collapse to one row per locus: repNames ordered by overlap bp, descending
awk -F'\t' -v OFS='\t' '
  { key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7
    if (!(key in seen)) { seen[key]=1; order[++n]=key }
    if ($14 > 0) { ov[key"\t"$11] += $14; cls[key"\t"$11]=$12
                   if ($6 == $13) conc[key"\t"$11]=1 }
  }
  END{
    for (i=1; i<=n; i++) {
      k=order[i]; m=0; delete nm; delete bp
      for (j in ov) { split(j,p,"\t"); kk=p[1]"\t"p[2]"\t"p[3]"\t"p[4]"\t"p[5]"\t"p[6]"\t"p[7]
                      if (kk==k) { nm[++m]=p[8]; bp[m]=ov[j] } }
      # insertion sort on overlap bp, descending
      for (a=2; a<=m; a++) { vn=nm[a]; vb=bp[a]; b=a-1
        while (b>=1 && bp[b]<vb) { nm[b+1]=nm[b]; bp[b+1]=bp[b]; b-- }
        nm[b+1]=vn; bp[b+1]=vb }
      s=""; for (a=1; a<=m; a++) s = s (a>1?"|":"") nm[a]
      print k, (m?s:"none")
    }
  }' "$W/annot_raw.tsv" | sort -k1,1 -k2,2n > "$W/union.bed8"

cut -f1-6,8 "$W/union.bed8" > "$OUT/AluACA_union_nr.bed"
{ printf "chrom\tstart\tend\tname\tscore\tstrand\tsource\trepeat_family\n"; cat "$W/union.bed8"; } \
  > "$OUT/AluACA_union_nr.tsv"

echo "[09] summary"
awk -F'\t' 'NR>1{
    n++; if ($8=="none") none++; else if ($8 ~ /^Alu/) alu++; else other++
    split($8,a,"|"); if ($8!="none") fam[a[1]]++
    src[$7"\t"($8 ~ /^Alu/ ? "Alu" : ($8=="none"?"none":"other"))]++ }
  END{ printf "  loci: %d   Alu-dominant: %d   other repeat: %d   no repeat: %d\n",n,alu,other,none
       print "  by source:"
       for (k in src) { split(k,p,"\t"); printf "    %-14s %-6s %d\n",p[1],p[2],src[k] } }' \
  "$OUT/AluACA_union_nr.tsv"
echo "  top repNames:"
awk -F'\t' 'NR>1 && $8!="none"{split($8,a,"|"); print a[1]}' "$OUT/AluACA_union_nr.tsv" \
  | sort | uniq -c | sort -rn | head -12 | sed 's/^/    /'
echo "[09] rewrote $OUT/AluACA_union_nr.bed (BED6+1) and $OUT/AluACA_union_nr.tsv"
