#!/usr/bin/env bash
# Step 01 - Fetch the AluACA RNA gene sequences from NCBI as FASTA.
#
# The full series is HE855917-HE856264 = 348 records = AluACA1-348, verified at
# both boundaries (HE855917 = AluACA1, HE856264 = AluACA348). HE855915/HE855916
# belong to an unrelated Cryptococcus laurentii submission, so the lower bound
# is a genuine edge. Starting at HE855919 instead silently drops AluACA1 and
# AluACA2 while still returning a plausible-looking 346.
#
# Gotcha: the square brackets in the Entrez [ACCN] field qualifier MUST be
# URL-encoded (%5B / %5D). Sending them raw returns an empty body with exit 0,
# which looks like a network failure but is not.
#
# efetch is retried: NCBI throttles bursts, and a throttled call returns an
# empty body rather than an HTTP error.

source "$(dirname "$0")/config.sh"

echo "[01] verifying accession range ..."
q() {
  curl -s --max-time 30 "$EUTILS/esearch.fcgi?db=nuccore&term=$1&retmode=json" \
    | grep -o '"count":"[0-9]*"' | head -1 | grep -o '[0-9]*'
}
echo "  $ACC_FROM              -> $(q "${ACC_FROM}%5BACCN%5D") record(s)"
echo "  $ACC_TO              -> $(q "${ACC_TO}%5BACCN%5D") record(s)"
range_n=$(q "${ACC_FROM}%3A${ACC_TO}%5BACCN%5D")
echo "  ${ACC_FROM}:${ACC_TO} -> ${range_n} record(s)"

expect=$(( $(echo "$ACC_TO" | tr -d 'A-Z') - $(echo "$ACC_FROM" | tr -d 'A-Z') + 1 ))

echo "[01] fetching FASTA ..."
ok=0
for attempt in 1 2 3; do
  r=$(curl -s --max-time 60 \
    "$EUTILS/esearch.fcgi?db=nuccore&term=${ACC_FROM}%3A${ACC_TO}%5BACCN%5D&usehistory=y&retmax=500&tool=$NCBI_TOOL&retmode=json") || true
  wenv=$(echo "$r" | grep -o '"webenv":"[^"]*"'   | sed 's/.*:"//;s/"//' || true)
  qk=$(  echo "$r" | grep -o '"querykey":"[^"]*"' | sed 's/.*:"//;s/"//' || true)
  if [ -z "$wenv" ] || [ -z "$qk" ]; then
    echo "  attempt $attempt: esearch returned no history handle; retrying" >&2
    sleep 5; continue
  fi
  curl -s --max-time 180 \
    "$EUTILS/efetch.fcgi?db=nuccore&WebEnv=$wenv&query_key=$qk&rettype=fasta&retmode=text&retmax=500&tool=$NCBI_TOOL" \
    > "$FASTA" || true
  n=$(grep -c '^>' "$FASTA" 2>/dev/null || true); n=${n:-0}
  if [ "$n" -eq "$expect" ]; then ok=1; break; fi
  echo "  attempt $attempt: got $n records, expected $expect; retrying" >&2
  sleep 5
done
[ "$ok" -eq 1 ] || { echo "[01] FAILED: could not fetch $expect records" >&2; exit 1; }

# --- integrity checks ---------------------------------------------------
echo "[01] records=$n  expected=$expect"
grep '^>' "$FASTA" | sed 's/^>\([A-Z0-9]*\)\..*/\1/' | tr -d 'A-Z' | sort > "$WORK/got.acc"
seq "$(echo "$ACC_FROM" | tr -d 'A-Z')" "$(echo "$ACC_TO" | tr -d 'A-Z')" | sort > "$WORK/want.acc"
echo "  missing accessions: $(comm -13 "$WORK/got.acc" "$WORK/want.acc" | wc -l)"
echo "  empty sequences:    $(awk '/^>/{if(n==0&&seen)e++;n=0;seen=1;next}{n+=length($0)}END{if(n==0)e++;print e+0}' "$FASTA")"
echo "  AluACA ID range:    $(grep '^>' "$FASTA" | sed -E 's/.*AluACA([0-9]+).*/\1/' | sort -n | sed -n '1p;$p' | paste -sd-)"
rm -f "$WORK/got.acc" "$WORK/want.acc"
echo "[01] wrote $FASTA"
