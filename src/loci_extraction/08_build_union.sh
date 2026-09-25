#!/usr/bin/env bash
# Step 08 - Non-redundant union of the two AluACA locus sets.
#
#   set A: the 344 placed Jady AluACAs           (AluACA_hg38.bed)
#   set B: the 543 NapRNAdb Alu/L1 ACA loci      (napRNA_Alu_L1_ACA.csv)
#
# 121 loci are shared. Because the deposited AluACA sequences are 3'-partial
# (median 79 nt) while the CSV carries the ~160 nt full element, the two sets
# describe the same locus at different extents. Where they overlap the union
# keeps the LONGER interval and joins both identifiers, so no locus is
# represented by a truncated interval.
#
# Two optional plausibility filters on the NapRNAdb side, both OFF by default
# so a rerun reproduces the committed catalogue. Whatever they leave behind,
# the summary at the end warns about intervals longer than a full-length Alu.
#
#   --max-len N   drop intervals >= N nt. The NapRNAdb novel-ACA entries run up
#                 to 6,859 nt -- a host intron or a LINE rather than an RNA, and
#                 the longest span 3-8 separate RepeatMasker elements while being
#                 annotated after just one of them. An H/ACA RNA is ~100-200 nt
#                 and a full-length Alu ~300, so --max-len 300 removes 53 CSV
#                 rows; --max-len 1000 removes only the 9 worst.
#                 Applied AFTER the intersect, because both sets have lengths:
#                 at a shared locus the kb-long partner is discarded rather than
#                 the locus, so no AluACA is lost either way.
#
#   --min-cov N   drop CSV rows whose Coverage column is below N. Applied BEFORE
#                 the intersect, because coverage is a property of the NapRNAdb
#                 call alone -- the Jady deposits have none -- so a row failing
#                 its own support threshold is not a locus to fall back to. A
#                 shared locus survives as jady_aluaca via the -v intersects
#                 below, so again no AluACA is lost. A missing or non-numeric
#                 Coverage counts as 0 and is dropped whenever the filter is on.
#
# Why this matters downstream: step 08 writes one FASTA record per interval, and
# the chimeric pipeline picks a read's guide by best bowtie2 hit across the whole
# catalogue. A record's chance of winning scales with how much sequence it holds,
# so a multi-kb interval collects guide assignments by area rather than identity.
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

# CSV columns: 1 Browser, 2 napRNA ID, 3 Chrom, 4 Start, 5 End, 6 Coverage,
# 7 Strand, 8 Length. The BED score column carries NapRNAdb's Length, not the
# coverage, which is consumed here.
awk -F',' -v OFS='\t' -v MINCOV="$MINCOV" -v N="$W/lowcov.n" '
  NR>1 { if (MINCOV > 0 && $6+0 < MINCOV) { dropped++; next }
         print $3,$4,$5,$2,$8,$7 }
  END  { print dropped+0 > N }' "$CSV" | sort -k1,1 -k2,2n > "$W/csv.bed"
sort -k1,1 -k2,2n "$OUT/AluACA_hg38.bed" > "$W/aluaca_raw.bed"

# collapse identical-interval AluACA duplicates, joining the IDs with "|"
"$BEDTOOLS" merge -i "$W/aluaca_raw.bed" -s -c 4,6 -o distinct,distinct \
  | awk -v OFS='\t' '{gsub(",","|",$4); print $1,$2,$3,$4,"0",$5}' \
  | sort -k1,1 -k2,2n > "$W/aluaca.bed"

# shared loci: pair them up, then keep whichever interval is longer
"$BEDTOOLS" intersect -a "$W/aluaca.bed" -b "$W/csv.bed" -wa -wb -s > "$W/pairs.tsv"
# ... except that with --max-len a kb-long partner is never the one kept: fall
# back to the shorter interval instead of dropping a locus with a real AluACA.
awk -F'\t' -v OFS='\t' -v MAX="$MAXLEN" '
  { la=$3-$2; lc=$9-$8
    jady = (la > lc)
    if (MAX > 0) {
      if ( jady && la >= MAX && lc < MAX) jady = 0
      if (!jady && lc >= MAX && la < MAX) jady = 1
      if (la >= MAX && lc >= MAX) next
    }
    if (jady) print $1,$2,$3,$4"|"$10,0,$6,"both_jady"
    else      print $7,$8,$9,$4"|"$10,0,$12,"both_naprnadb" }' \
  "$W/pairs.tsv" > "$W/shared.bed7"

# the rest of each set, untouched
"$BEDTOOLS" intersect -a "$W/aluaca.bed" -b "$W/csv.bed" -v -s \
  | awk -v OFS='\t' -v MAX="$MAXLEN" 'MAX<=0||($3-$2)<MAX{print $1,$2,$3,$4,0,$6,"jady_aluaca"}'   > "$W/jady_only.bed7"
"$BEDTOOLS" intersect -a "$W/csv.bed" -b "$W/aluaca.bed" -v -s \
  | awk -v OFS='\t' -v MAX="$MAXLEN" 'MAX<=0||($3-$2)<MAX{print $1,$2,$3,$4,0,$6,"naprnadb_only"}' > "$W/csv_only.bed7"

cat "$W/shared.bed7" "$W/jady_only.bed7" "$W/csv_only.bed7" \
  | sort -k1,1 -k2,2n > "$W/union.bed7"

cut -f1-6 "$W/union.bed7" > "$OUT/AluACA_union_nr.bed"
{ printf "chrom\tstart\tend\tname\tscore\tstrand\tsource\n"; cat "$W/union.bed7"; } \
  > "$OUT/AluACA_union_nr.tsv"

# stranded sequence straight from hg38, so both sources are on equal footing.
# Headers follow the snoRNA.txt.fa convention -- ">name.idN", one unwrapped
# uppercase sequence line -- so the file concatenates with that catalogue.
# Numbering starts above the highest id there (2089) so the two never collide,
# and "|" in a joined name becomes "_" because the format is a single token.
"$BEDTOOLS" getfasta -fi "$HG38_FA" -bed "$OUT/AluACA_union_nr.bed" -s -name \
  | awk -v BASE="$FASTA_ID_BASE" '
      /^>/ { if (n) print h "\n" s
             h = $0; sub(/::.*$/, "", h); gsub(/\|/, "_", h)
             h = h ".id" (BASE + ++i); s = ""; n = 1; next }
      { s = s toupper($0) }
      END { if (n) print h "\n" s }' > "$OUT/AluACA_union_nr.fasta"

# NB: config.sh sets -e, so these stay as explicit ifs -- a trailing
# `[ test ] && cmd` would make this block exit 1 whenever the test is false.
_filters=()
if [ "${MAXLEN:-0}" -gt 0 ]; then
  _filters+=("intervals >= $MAXLEN nt dropped")
fi
if [ "${MINCOV:-0}" -gt 0 ]; then
  _filters+=("NapRNAdb coverage < $MINCOV dropped: $(cat "$W/lowcov.n") row(s)")
fi
if [ "${#_filters[@]}" -gt 0 ]; then
  # ${arr[*]} joins on the FIRST character of IFS only, so build the "; "
  # separator explicitly rather than getting "a;b".
  _msg=$(printf '%s; ' "${_filters[@]}"); _msg=${_msg%'; '}
  echo "[08] union summary  ($_msg)"
else
  echo "[08] union summary  (no filters; --max-len / --min-cov are both 0)"
fi
cut -f7 "$W/union.bed7" | sort | uniq -c | awk '{printf "  %-16s %s\n",$2,$1}'
echo "  union intervals:               $(wc -l < "$OUT/AluACA_union_nr.bed")"
echo "  fasta records:                 $(grep -c '^>' "$OUT/AluACA_union_nr.fasta")"
echo "  residual same-strand overlaps: $(( $(wc -l < "$OUT/AluACA_union_nr.bed") - $("$BEDTOOLS" merge -s -i "$OUT/AluACA_union_nr.bed" | wc -l) ))"
echo "  duplicate names:               $(( $(cut -f4 "$OUT/AluACA_union_nr.bed" | wc -l) - $(cut -f4 "$OUT/AluACA_union_nr.bed" | sort -u | wc -l) ))"
echo "  AluACA sequence not covered by the chosen interval:"
awk -F'\t' '($9-$8)>=($3-$2) && ($2<$8 || $3>$9) {printf "    %-22s aluaca %s:%d-%d  kept %s:%d-%d\n",$4,$1,$2,$3,$7,$8,$9}' "$W/pairs.tsv"

# An H/ACA RNA is ~100-200 nt and a full-length Alu ~300. Anything past that
# becomes a FASTA record that wins guide assignment on length alone, so say so
# rather than letting it through quietly.
_over=$(awk -F'\t' '($3-$2)>=300' "$W/union.bed7" | wc -l)
if [ "$_over" -gt 0 ]; then
  echo "  !! $_over interval(s) >= 300 nt, longer than a full-length Alu."
  echo "     Each becomes one FASTA record and collects guide assignments by"
  echo "     sequence length rather than identity.  Re-run with --max-len 300."
  awk -F'\t' '($3-$2)>=300{print $3-$2"\t"$4"\t"$1":"$2"-"$3"\t"$7}' "$W/union.bed7" \
    | sort -k1,1nr | head -5 \
    | awk -F'\t' '{printf "       %6d nt  %-30s %s  (%s)\n",$1,$2,$3,$4}'
  if [ "$_over" -gt 5 ]; then
    echo "       ... and $((_over - 5)) more"
  fi
fi
