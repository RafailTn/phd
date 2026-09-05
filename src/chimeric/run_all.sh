#!/usr/bin/env bash
# Drive the chimeric eCLIP pipeline end to end, for one or more arms.
#
#   bash src/chimeric/run_all.sh                  # the default arm, both samples
#   bash src/chimeric/run_all.sh arm1 arm3        # two arms, reusing the index
#   bash src/chimeric/run_all.sh arm0d            # the workstation run again, dense
#   bash src/chimeric/run_all.sh --species hg19 --source merged SRR30692552
#   bash src/chimeric/run_all.sh arm3 --only annotate,compare
#
# Steps, in order. Each is skipped with a notice when its inputs are already
# there, so re-running is cheap and safe:
#
#   fastq     fetch_fastq.sh     download the runs into $WORK
#   refs      fetch_refs.sh      download genome + GENCODE + rmsk into $REF
#   index     build_indices.sh   build the STAR indices
#   run       run_chimeras.sh    the pipeline itself, per sample
#   annotate  annotate_chimeras.py   add gene/repeat context
#   compare   compare_to_published.py  against the published hg19 chimeras
#   report    make_report.py     regenerate RESULTS.md
#
#   --only a,b    run only these steps
#   --skip a,b    run everything except these
#
# The reproduction matrix these arms make up:
#
#            | plain snoRNA.txt.fa      | merged AluACA+snoRNA
#   ---------+--------------------------+--------------------------
#   hg38     | arm2                     | arm0 (sparse) / arm0d (dense)
#   hg19     | arm1  <- reproduction    | arm3  <- build-stability
#
# arm1 is the number to judge against the published 45,810 genomic chimeras in
# DKC1_IP.snoRNA.hg19.chimeras.csv. arm3 says whether the AluACA findings survive
# a change of build. arm0d re-runs the workstation's own configuration densely
# and so measures the sparse-vs-dense term every sparse arm carries.
#
# Every path is settable by flag or environment variable; run --help for the
# full list.
set -euo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Bare armN tokens are pulled out before config.sh is sourced, because each arm
# needs its own species, catalogue and index density: config.sh derives those
# with ${VAR:-default}, so a second source in this process would keep the first
# arm's values. Re-running this script per arm gives each a clean environment.
ARMS=() ; REST=() ; _prev=
for _a in "$@"; do
  case "$_prev,$_a" in
    # the value of --arm is config.sh's business, not a token to dispatch on
    '--arm,'*)   REST+=("$_a") ;;
    *,arm[0-9]*) ARMS+=("$_a") ;;
    *)           REST+=("$_a") ;;
  esac
  _prev=$_a
done
if [ "${#ARMS[@]}" -gt 0 ]; then
  for _arm in "${ARMS[@]}"; do
    bash "$SRC/run_all.sh" --arm "$_arm" ${REST[@]+"${REST[@]}"}
  done
  echo ALL_ARMS_COMPLETE
  exit 0
fi

source "$SRC/config.sh"

# --- run_all's own flags, from the arguments config.sh passed through --------
ONLY= ; SKIP= ; SAMPLES=()
while [ "$#" -gt 0 ]; do
  case "$1" in
    --only) ONLY="$2"; shift 2 ;;
    --skip) SKIP="$2"; shift 2 ;;
    -*)     echo "unknown flag: $1  (see --help)" >&2; exit 1 ;;
    *)      SAMPLES+=("$1"); shift ;;
  esac
done
# No sample named: the IP and its input control.
[ "${#SAMPLES[@]}" -eq 0 ] && read -r -a SAMPLES <<<"$SRRS"

want () {  # want <step> -- is this step in scope?
  [ -n "$ONLY" ] && { case ",$ONLY," in *",$1,"*) return 0 ;; *) return 1 ;; esac; }
  [ -n "$SKIP" ] && { case ",$SKIP," in *",$1,"*) return 1 ;; esac; }
  return 0
}

say () { echo "######## $* $(date '+%F %T') ########"; }

# --- steps ------------------------------------------------------------------
# Each reads the configuration for the arm currently being set up, so the
# per-arm re-source below is what switches species, catalogue and density.

step_fastq () {
  want fastq || return 0
  local missing=0
  for srr in "${SAMPLES[@]}"; do [ -s "$WORK/$srr.fastq.gz" ] || missing=1; done
  if [ "$missing" -eq 0 ]; then echo "[fastq] all present in $WORK"; return 0; fi
  say "fetch fastq"
  bash "$SRC/fetch_fastq.sh" "${SAMPLES[@]}"
}

step_refs () {
  want refs || return 0
  if [ -f "$GENOME_FA" ] && [ -f "$GENCODE" ] && [ -s "$RMSK_BED" ]; then
    echo "[refs] $SPECIES references present"; return 0
  fi
  say "fetch $SPECIES references"
  bash "$SRC/fetch_refs.sh" "$SPECIES"
}

step_index () {
  want index || return 0
  say "build $SPECIES indices (sparseD=$SPARSE_D)"
  bash "$SRC/build_indices.sh"
}

step_run () {
  want run || return 0
  for srr in "${SAMPLES[@]}"; do
    say "$ARM $srr run"
    bash "$SRC/run_chimeras.sh" "$srr" --outdir "$OUT/$ARM/$srr"
  done
}

step_annotate () {
  want annotate || return 0
  # The BEDs are optional context: without them the flag columns are simply
  # absent, which is better than refusing to annotate at all.
  local bed_args=()
  [ -s "$GUIDE_BED" ] && bed_args+=(--source-bed "$GUIDE_BED") \
    || echo "[annotate] no guide-locus BED at $GUIDE_BED, skipping that flag"
  [ -s "$RMSK_BED" ] && bed_args+=(--rmsk "$RMSK_BED") \
    || echo "[annotate] no RepeatMasker BED at $RMSK_BED, skipping that flag"
  for srr in "${SAMPLES[@]}"; do
    say "$ARM $srr annotate"
    "$PYTHON" "$SRC/annotate_chimeras.py" \
      --outdir "$OUT/$ARM/$srr" --uid "$srr" \
      --gtag "$SPECIES" --tags "$(echo "$TARGET_TAGS" | tr ' ' ','),$SPECIES" \
      --gtf "$GENCODE" --alu-fasta "$ALU_FASTA" \
      --bedtools "$BEDTOOLS" \
      ${bed_args[@]+"${bed_args[@]}"} \
      --out "$OUT/$ARM/$srr.annotated.tsv"
  done
}

step_compare () {
  want compare || return 0
  if [ ! -s "$PUBLISHED" ]; then
    echo "[compare] skipped: no published chimeras CSV at $PUBLISHED  (set --published)"
    return 0
  fi
  say "$ARM compare vs published"
  # A label describing what was actually run, not what the arm table says: the
  # density can be overridden with --sparse / --dense.
  local density=sparse; [ "$SPARSE_D" = 1 ] && density=dense
  "$PYTHON" "$SRC/compare_to_published.py" \
    --outdir "$OUT/$ARM/$IP" --uid "$IP" --gtag "$SPECIES" \
    --published "$PUBLISHED" \
    --label "$ARM ($SPECIES + $SOURCE, $density)" \
    --out "$OUT/$ARM.vs_published.txt"
}

step_report () {
  want report || return 0
  # The report contrasts the IP against its input control, so it needs both.
  if [ "${#SAMPLES[@]}" -lt 2 ]; then
    echo "[report] skipped: needs the IP and its input control, got ${SAMPLES[*]}"
    return 0
  fi
  say "$ARM report"
  "$PYTHON" "$SRC/make_report.py" \
    --arm "$ARM" --resdir "$OUT/$ARM" \
    --ip "${SAMPLES[0]}" --input "${SAMPLES[1]}" \
    --published "$PUBLISHED" --rmsk "$RMSK_BED" \
    --bedtools "$BEDTOOLS" --gtag "$SPECIES" \
    --out "$OUT/RESULTS.md"
}

# --- drive ------------------------------------------------------------------
step_fastq; step_refs; step_index
step_run;   step_annotate; step_compare; step_report
echo "######## $ARM complete $(date '+%F %T') ########"
