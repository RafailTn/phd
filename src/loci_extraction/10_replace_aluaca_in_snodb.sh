#!/usr/bin/env bash
# Step 10 - snoDB catalogue with its AluACA tier swapped for our union set.
#
# Takes the full snoDB 2.0 TSV export, drops every row whose box_type is
# "AluACA" (352 of them - fixed-width windows that are largely the Jady set
# re-expressed), and substitutes the 765-locus union built by step 08.
#
# COORDINATES: the snoDB export is 1-based inclusive (verified: end-start+1
# equals its own `length` column for all 2123 rows). The union is BED, 0-based
# half-open. The output is BED convention throughout - snoDB starts are
# decremented on the way in - so columns 1-6 are a valid BED6 and bedtools can
# read it directly.
#
# Columns: chrom start end name score strand origin box_type
#   origin    snoDB | AluACA_union
#   box_type  snoDB's own value, or "Alu-ACA" for every union locus
#
# A few snoDB rows that escaped that filter still describe union loci: they are
# classified "H/ACA" rather than "AluACA" but sit on the same interval and
# strand. Any retained row that overlaps a union locus same-strand is dropped
# too, so no locus is counted twice.
#
# Names must be unique for the file to work as a reference: snoDB gene_name is
# empty for 733 rows and duplicated for 50, so a row falls back to its
# snodb_id, and any remaining collision gets the snodb_id appended.
#
# Writes snoDB_with_AluACA_union.tsv and .fasta to $OUT. Sequence is pulled
# stranded from hg38 for every row, so both origins are on equal footing; the
# snoDB rows are then checked against snoDB's own `sequence` column, which is
# also what validates the 1-based -> 0-based conversion above.

source "$(dirname "$0")/config.sh"

UNION="$OUT/AluACA_union_nr.bed"
DEST="$OUT/snoDB_with_AluACA_union.tsv"
FA="$OUT/snoDB_with_AluACA_union.fasta"

for f in "$SNODB_TSV" "$UNION"; do
  [ -f "$f" ] || { echo "MISSING: $f" >&2; exit 1; }
done

W="$WORK/snodb_merge"; mkdir -p "$W"

# --- snoDB rows to keep, converted to BED coordinates -------------------
# Pass 1 counts each gene_name so pass 2 knows which ones collide.
awk -F'\t' -v OFS='\t' '
  NR==1 { next }
  $19 == "AluACA" { dropped++; next }
  $13 == "" || $14 == "" { noloc++; next }
  { rows[++n] = $0; if ($17 != "") freq[$17]++ }
  END {
    for (i = 1; i <= n; i++) {
      split(rows[i], c, "\t")
      name = (c[17] != "" && freq[c[17]] == 1) ? c[17] : \
             (c[17] != "" ? c[17] "_" c[1] : c[1])
      print c[13], c[14] - 1, c[15], name, 0, c[16], "snoDB", c[19]
      if (c[23] != "") print name "\t" toupper(c[23]) > SEQREF
    }
    printf "  snoDB rows kept:        %d\n", n              > "/dev/stderr"
    printf "  AluACA rows dropped:    %d\n", dropped        > "/dev/stderr"
    if (noloc) printf "  skipped, no coords:     %d\n", noloc > "/dev/stderr"
  }' SEQREF="$W/snodb_seq.tsv" "$SNODB_TSV" > "$W/keep.bed8"

# --- the union, relabelled ----------------------------------------------
awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,0,$6,"AluACA_union","Alu-ACA"}' \
  "$UNION" > "$W/union.bed8"

# --- drop retained rows that duplicate a union locus ---------------------
sort -k1,1 -k2,2n "$W/keep.bed8" > "$W/keep.sorted.bed8"
"$BEDTOOLS" intersect -a "$W/keep.sorted.bed8" -b "$W/union.bed8" -wa -wb -s \
  > "$W/redundant.tsv"
"$BEDTOOLS" intersect -a "$W/keep.sorted.bed8" -b "$W/union.bed8" -v -s \
  > "$W/keep.final.bed8"

if [ -s "$W/redundant.tsv" ]; then
  echo "  dropped as duplicates of a union locus:"
  awk -F'\t' '{printf "    %-22s %s:%s-%s(%s) %-8s <-> %s\n",$4,$1,$2,$3,$6,$8,$12}' \
    "$W/redundant.tsv"
fi

{ printf "chrom\tstart\tend\tname\tscore\tstrand\torigin\tbox_type\n"
  sort -k1,1 -k2,2n "$W/keep.final.bed8" "$W/union.bed8"
} > "$DEST"

echo "[10] $DEST"
awk -F'\t' 'NR>1{o[$7]++; b[$8]++} END{
    print "  by origin:";   for (k in o) printf "    %-14s %d\n", k, o[k]
    print "  by box_type:"; for (k in b) printf "    %-14s %d\n", (k==""?"(empty)":k), b[k] }' "$DEST"
echo "  rows:             $(( $(wc -l < "$DEST") - 1 ))"
echo "  duplicate names:  $(( $(tail -n +2 "$DEST" | cut -f4 | wc -l) - $(tail -n +2 "$DEST" | cut -f4 | sort -u | wc -l) ))"

# --- stranded sequence from hg38 ----------------------------------------
tail -n +2 "$DEST" | cut -f1-6 > "$W/all.bed"
"$BEDTOOLS" getfasta -fi "$HG38_FA" -bed "$W/all.bed" -s -name \
  | sed -E 's/^>([^:]+)::/>\1 /' > "$FA"
echo "[10] $FA"
echo "  records:          $(grep -c '^>' "$FA")"

# snoDB ships its own sequence for most rows; agreement confirms both the
# coordinates and the 1-based -> 0-based conversion.
awk '/^>/{n=substr($1,2); next}{s[n]=s[n]$0} END{for(k in s) print k"\t"toupper(s[k])}' "$FA" \
  > "$W/extracted_seq.tsv"
awk -F'\t' '
  NR==FNR { ref[$1]=$2; next }
  ($1 in ref) { n++; if ($2 == ref[$1]) same++; else { diff++; if (diff<=5) bad = bad "    " $1 "\n" } }
  END { printf "  vs snoDB `sequence` column: %d compared, %d identical, %d differ\n", n, same, diff
        if (diff) printf "%s", bad }' \
  "$W/snodb_seq.tsv" "$W/extracted_seq.tsv" 
