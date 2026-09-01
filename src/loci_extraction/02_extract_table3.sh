#!/usr/bin/env bash
# Step 02 - Extract "Supplemental Table 3. List of the host genes of human
# intron-encoded AluACA RNAs" from the Jady supplemental PDF.
#
# Three parsing traps in this PDF, each of which silently drops entries:
#   1. Page-break form feeds (\f) prefix the first line of each page, so an
#      anchored /^AluACA/ match misses those lines (AluACA71,148,224,301).
#   2. The AluACA105 line has a stray trailing "AluACA106," from the original
#      Word layout, and the real 106 entry is typo'd as "ALuACA106" (capital L).
#   3. Entries are "AluACA<N>, <SYMBOL>, <description>" but some omit the
#      symbol and give only a description.
# Handling 1+2 recovers the full 348.

source "$(dirname "$0")/config.sh"

echo "[02] pdftotext ..."
pdftotext -layout "$PDF" "$WORK/pdf.txt"

# Table 3 begins at the line matching its caption; take everything after it.
start=$(grep -n "Supplemental Table 3" "$WORK/pdf.txt" | head -1 | cut -d: -f1)
[ -n "$start" ] || { echo "Table 3 caption not found" >&2; exit 1; }
echo "[02] Table 3 caption at line $start"

awk -v s="$start" 'NR>s' "$WORK/pdf.txt" \
  | tr -d '\f' \
  | grep -iE "^AluACA[0-9]+," \
  | sed -E 's/[[:space:]]*AluACA106,[[:space:]]*$//I' \
  | sed -E 's/^ALuACA/AluACA/I' \
  > "$WORK/table3.txt"

n=$(wc -l < "$WORK/table3.txt")
echo "[02] entries: $n"
awk -F, '{id=$1; sub(/[Aa][Ll][Uu][Aa][Cc][Aa]/,"",id); have[id]=1}
 END{m=0; for(i=1;i<=348;i++) if(!(i in have)){printf "  MISSING AluACA%s\n",i; m++}
     print "  missing total: " m}' "$WORK/table3.txt"
echo "  duplicate IDs: $(cut -d, -f1 "$WORK/table3.txt" | sort | uniq -d | tr '\n' ' ')"
echo "  'No apparent host gene': $(grep -ci 'No apparent host gene' "$WORK/table3.txt")"
echo "[02] wrote $WORK/table3.txt"
