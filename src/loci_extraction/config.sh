#!/usr/bin/env bash
# Shared configuration for the AluACA -> hg38 mapping pipeline.
# Sourced by every step:  source "$(dirname "$0")/config.sh"
#
# Every path is settable three ways, in increasing precedence:
#   1. built-in default   (project = the directory containing scripts/)
#   2. environment variable
#   3. command-line flag
#
# Because this file is sourced, the flags below work on any step script as
# well as on run_all.sh:
#
#   bash scripts/run_all.sh --proj /data/aluaca --hg38-dir /ref/hg38 \
#                           --rmsk /ref/rmsk/hg38_rmsk.tsv.gz --out /results
#   bash scripts/09_add_repeat_family.sh --rmsk /ref/rmsk/hg38_rmsk.tsv.gz
#
# Run any script with --help for the full list.

set -euo pipefail

_cfg_usage() {
  cat <<'USAGE'
AluACA -> hg38 pipeline - configuration flags

  --proj DIR         project directory (inputs + default output location)
                     default: the parent of scripts/
  --out DIR          where deliverables are written        [default: $PROJ]
  --work DIR         scratch directory for intermediates   [default: $OUT/work_map]

  --hg38-dir DIR     directory holding the genome + GENCODE
  --hg38-fa FILE     genome FASTA (must be samtools-indexed, .fai alongside)
  --gencode FILE     GENCODE annotation GTF(.gz)
  --rmsk FILE        UCSC RepeatMasker table (tab-separated, gzipped ok)
                     required by step 09 only
  --max-len N        drop union intervals of N nt or longer
                     default 0 = no length filter

  --csv FILE         napRNAdb Alu/L1 ACA CSV
  --pdf FILE         Jady et al. supplemental PDF
  --fasta FILE       deposited-sequence FASTA (written by step 01)

  --bin DIR          directory containing bedtools + python3
                     default: deps/.pixi/envs/default/bin if present, else PATH
  --bedtools PATH    bedtools executable
  --python PATH      python3 executable

  --acc-from ACC     first EMBL accession            [default: HE855917]
  --acc-to ACC       last EMBL accession             [default: HE856264]

  -h, --help         show this and exit

Environment variables of the same name in upper snake case (PROJ, OUT, WORK,
HG38_DIR, HG38_FA, GENCODE, RMSK, MAXLEN, CSV, PDF, FASTA, BIN, BEDTOOLS,
PYTHON, ACC_FROM, ACC_TO) are honoured as defaults; flags win over them.
USAGE
}

# --- flag parsing -------------------------------------------------------
# Unrecognised arguments are left alone so individual steps can define their
# own; only the shared ones are consumed here.
_cfg_rest=()
while [ "$#" -gt 0 ]; do
  case "$1" in
    --proj)      PROJ="$2";      shift 2 ;;
    --out)       OUT="$2";       shift 2 ;;
    --work)      WORK="$2";      shift 2 ;;
    --hg38-dir)  HG38_DIR="$2";  shift 2 ;;
    --hg38-fa)   HG38_FA="$2";   shift 2 ;;
    --gencode)   GENCODE="$2";   shift 2 ;;
    --rmsk)      RMSK="$2";      shift 2 ;;
    --max-len)   MAXLEN="$2";    shift 2 ;;
    --csv)       CSV="$2";       shift 2 ;;
    --pdf)       PDF="$2";       shift 2 ;;
    --fasta)     FASTA="$2";     shift 2 ;;
    --bin)       BIN="$2";       shift 2 ;;
    --bedtools)  BEDTOOLS="$2";  shift 2 ;;
    --python)    PYTHON="$2";    shift 2 ;;
    --acc-from)  ACC_FROM="$2";  shift 2 ;;
    --acc-to)    ACC_TO="$2";    shift 2 ;;
    -h|--help)   _cfg_usage; exit 0 ;;
    *)           _cfg_rest+=("$1"); shift ;;
  esac
done
set -- ${_cfg_rest[@]+"${_cfg_rest[@]}"}

# --- project layout -----------------------------------------------------
# The project root is found by walking up from this file until a directory
# holding the pipeline's inputs turns up, so the scripts can be nested at any
# depth (scripts/, scripts/loci_extraction/, ...) without breaking.
_cfg_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -z "${PROJ:-}" ]; then
  PROJ="$_cfg_dir"
  while [ "$PROJ" != "/" ]; do
    if [ -e "$PROJ/napRNA_Alu_L1_ACA.csv" ] || [ -d "$PROJ/deps" ]; then break; fi
    PROJ="$(dirname "$PROJ")"
  done
  # nothing recognisable found: fall back to the parent of the scripts dir
  [ "$PROJ" = "/" ] && PROJ="$(dirname "$_cfg_dir")"
fi
OUT="${OUT:-$PROJ}"
WORK="${WORK:-$OUT/work_map}"

# --- external reference data -------------------------------------------
# Filenames are guessed inside HG38_DIR when not given explicitly, so a server
# carrying a different GENCODE release or assembly filename still works.
HG38_DIR="${HG38_DIR:-$HOME/Downloads/hg38}"
if [ -z "${HG38_FA:-}" ]; then
  HG38_FA="$HG38_DIR/GRCh38.primary_assembly.genome.fa"
  [ -e "$HG38_FA" ] || HG38_FA="$(ls -1 "$HG38_DIR"/*.{fa,fasta} 2>/dev/null | head -1 || true)"
fi
if [ -z "${GENCODE:-}" ]; then
  GENCODE="$HG38_DIR/gencode.v47.primary_assembly.annotation.gtf.gz"
  [ -e "$GENCODE" ] || GENCODE="$(ls -1 "$HG38_DIR"/gencode*.gtf.gz 2>/dev/null | head -1 || true)"
fi
RMSK="${RMSK:-$HOME/Downloads/transposon_proj/data/hg38_rmsk.gtf.gz}"

# Optional length filter on the union (step 08): a handful of NapRNAdb entries
# span whole kilobases -- a host intron or a LINE, not an ACA RNA. 0 keeps
# everything; --max-len 1000 drops the 9 kb-long NapRNAdb-only loci.
MAXLEN="${MAXLEN:-0}"

# --- tools --------------------------------------------------------------
# Prefer the project's pixi env when it exists; otherwise fall back to PATH so
# the pipeline runs under a module system, conda env or plain container.
BIN="${BIN:-$PROJ/deps/.pixi/envs/default/bin}"
if [ -z "${BEDTOOLS:-}" ]; then
  if [ -x "$BIN/bedtools" ]; then BEDTOOLS="$BIN/bedtools"
  else BEDTOOLS="$(command -v bedtools || true)"; fi
fi
if [ -z "${PYTHON:-}" ]; then
  if [ -x "$BIN/python3" ]; then PYTHON="$BIN/python3"
  else PYTHON="$(command -v python3 || true)"; fi
fi

# --- inputs -------------------------------------------------------------
PDF="${PDF:-$PROJ/Supplemental_material.pdf}"      # Jady et al. supplemental tables
CSV="${CSV:-$PROJ/napRNA_Alu_L1_ACA.csv}"          # napRNAdb Alu/L1 ACA loci (hg38)
FASTA="${FASTA:-$OUT/AluACA_HE855917-HE856264.fasta}"   # produced by step 01

# --- NCBI ---------------------------------------------------------------
EUTILS="${EUTILS:-https://eutils.ncbi.nlm.nih.gov/entrez/eutils}"
ACC_FROM="${ACC_FROM:-HE855917}"
ACC_TO="${ACC_TO:-HE856264}"
# NCBI allows 3 requests/sec without an API key; steps that loop sleep 0.4s.
NCBI_TOOL="${NCBI_TOOL:-claude_code}"

export PROJ OUT WORK HG38_DIR HG38_FA GENCODE RMSK MAXLEN BIN BEDTOOLS PYTHON \
       PDF CSV FASTA EUTILS ACC_FROM ACC_TO NCBI_TOOL

mkdir -p "$WORK" "$OUT"

# Fail early with a clear message rather than midway through a step.
# RMSK is checked by step 09 itself, since only that step needs it.
_cfg_missing=0
_cfg_need() {   # _cfg_need <label> <path> <how-to-set>
  if [ -z "$2" ]; then
    echo "MISSING $1: could not be located automatically  ($3)" >&2; _cfg_missing=1
  elif [ ! -e "$2" ]; then
    echo "MISSING $1: $2  ($3)" >&2; _cfg_missing=1
  fi
}

[ -d "$HG38_DIR" ] || {
  echo "MISSING hg38 directory: $HG38_DIR  (set it with --hg38-dir)" >&2; _cfg_missing=1; }

_cfg_need "genome FASTA"  "$HG38_FA"  "--hg38-fa, or put a .fa in $HG38_DIR"
_cfg_need "GENCODE GTF"   "$GENCODE"  "--gencode, or put gencode*.gtf.gz in $HG38_DIR"
_cfg_need "bedtools"      "$BEDTOOLS" "--bedtools, --bin, or put it on PATH"
_cfg_need "python3"       "$PYTHON"   "--python, --bin, or put it on PATH"
[ -z "$HG38_FA" ] || _cfg_need "genome index" "${HG38_FA}.fai" \
  "run: samtools faidx '$HG38_FA'"

[ "$_cfg_missing" -eq 0 ] || {
  echo "" >&2
  echo "Every path is settable by flag or environment variable; see --help." >&2
  exit 1
}
