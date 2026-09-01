#!/usr/bin/env bash
# Run the whole AluACA -> hg38 pipeline end to end.
#
#   bash scripts/run_all.sh
#   bash scripts/run_all.sh --proj /data/aluaca --hg38-dir /ref/hg38 \
#                           --rmsk /ref/hg38_rmsk.tsv.gz --out /results
#
# Run with --help for every configurable path.
#
# Steps 01 and 04 need network access (NCBI E-utilities). Step 06 streams the
# 3.1 GB genome FASTA and is the slow one (~10-15 min); everything else is
# fast. Intermediates land in work_map/ and can be deleted afterwards.

source "$(dirname "$0")/config.sh"
D="$(dirname "$0")"

# config.sh has already parsed and exported the shared flags, so the steps are
# invoked with no arguments and inherit everything through the environment.
bash "$D/01_fetch_sequences.sh"
bash "$D/02_extract_table3.sh"
bash "$D/03_build_gencode_genes.sh"
"$PYTHON" "$D/04_resolve_symbols.py"
bash "$D/05_locate_in_genes.sh"
"$PYTHON" "$D/06_genome_search.py"
bash "$D/07_assemble_and_match.sh"
bash "$D/08_build_union.sh"
if [ -n "${RMSK:-}" ] && [ -f "$RMSK" ]; then
  bash "$D/09_add_repeat_family.sh"
else
  echo "[09] skipped: no RepeatMasker table at ${RMSK:-<unset>}  (set it with --rmsk)"
fi

echo
echo "done. deliverables:"
ls -la "$FASTA" \
       "$OUT/AluACA_hg38_coordinates.tsv" \
       "$OUT/AluACA_hg38.bed" \
       "$OUT/AluACA_unresolved.tsv" \
       "$OUT/AluACA_union_nr.bed" \
       "$OUT/AluACA_union_nr.tsv" \
       "$OUT/AluACA_union_nr.fasta"
echo
echo "to remove intermediates:  rm -rf $WORK"
