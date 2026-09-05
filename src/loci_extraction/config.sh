#!/usr/bin/env bash
# Shared configuration for the AluACA -> hg38 mapping pipeline.
# Sourced by every step:  source "$(dirname "$0")/config.sh"
#
# Every path is settable three ways, in increasing precedence:
#   1. built-in default   (project = the directory containing src/)
#   2. environment variable
#   3. command-line flag
#
# Because this file is sourced, the flags below work on any step script as
# well as on run_all.sh:
#
#   bash src/loci_extraction/run_all.sh --proj /data/aluaca --hg38-dir /ref/hg38 \
#                           --rmsk /ref/rmsk/hg38_rmsk.tsv.gz --out /results
#   bash src/loci_extraction/09_add_repeat_family.sh --rmsk /ref/rmsk/hg38_rmsk.tsv.gz
#
# Run any script with --help for the full list.

set -euo pipefail

_cfg_usage() {
  cat <<'USAGE'
AluACA -> hg38 pipeline - configuration flags

  --proj DIR         project directory (inputs + default output location)
                     default: the project root
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
  --polya-csv FILE   napRNAdb Alu/L1 polyA-pocket ACA CSV
                     required by csv_overlap_check.sh only
  --pdf FILE         Jady et al. supplemental PDF
  --fasta FILE       deposited-sequence FASTA (written by step 01)
  --snodb FILE       snoDB catalogue TSV
                     required by step 10 only
  --fasta-id-base N  first numeric id given to a union FASTA record in step 08
                     default 3000, above the highest id in snoRNA.txt.fa (2089)
                     so the two catalogues never collide

  --bin DIR          directory containing bedtools + python3
                     default: deps/.pixi/envs/default/bin if present, else PATH
  --bedtools PATH    bedtools executable
  --python PATH      python3 executable

  --acc-from ACC     first EMBL accession            [default: HE855917]
  --acc-to ACC       last EMBL accession             [default: HE856264]

  -h, --help         show this and exit

Environment variables of the same name in upper snake case (PROJ, OUT, WORK,
HG38_DIR, HG38_FA, GENCODE, RMSK, MAXLEN, CSV, POLYA_CSV, PDF, FASTA, SNODB_TSV,
FASTA_ID_BASE, BIN, BEDTOOLS, PYTHON, ACC_FROM, ACC_TO) are honoured as
defaults; flags win over them.
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
    --polya-csv) POLYA_CSV="$2";  shift 2 ;;
    --snodb)     SNODB_TSV="$2";  shift 2 ;;
    --fasta-id-base) FASTA_ID_BASE="$2"; shift 2 ;;
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
# depth (src/, src/loci_extraction/, ...) without breaking.
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

# --- input lookup -------------------------------------------------------
# An input is looked for at the project root, then in $PROJ/data (the more
# natural place to drop them on a shared server), then anywhere under $PROJ.
# The search never leaves the project: the current directory is deliberately
# not consulted, so running from somewhere else cannot silently pick up a
# same-named file. Flags override the lookup entirely.
_cfg_indexed=0
_cfg_index=
_cfg_build_index() {
  [ "$_cfg_indexed" = 1 ] && return 0
  _cfg_indexed=1
  # deps/ is a vendored conda environment of ~28k files and the STAR indices are
  # large; neither holds inputs, so pruning them keeps the walk to hundreds of
  # entries and stops a stray name in the env from matching. find does not
  # follow symlinks by default, so the walk cannot escape $PROJ either.
  _cfg_index=$(find "$PROJ" \
      \( -name .git -o -name deps -o -name __pycache__ -o -name '*_star_index' \) -prune \
      -o -type f -print 2>/dev/null)
}

_cfg_input() {   # _cfg_input <basename> <varname> [preferred dir ...]
  local base=$1 var=$2 d hit
  shift 2
  for d in "$@" "$PROJ" "$PROJ/data"; do
    if [ -e "$d/$base" ]; then printf -v "$var" '%s' "$d/$base"; return 0; fi
  done
  _cfg_build_index
  # Shallowest match wins, ties broken alphabetically, so which file is chosen
  # never depends on the order the filesystem happened to return.
  hit=$(printf '%s\n' "$_cfg_index" | awk -v b="$base" '
          { n = split($0, p, "/"); if (p[n] == b) print n "\t" $0 }' \
        | sort -k1,1n -k2,2 | head -1 | cut -f2-)
  # Report the root path when missing, so the error names where it should go.
  printf -v "$var" '%s' "${hit:-$PROJ/$base}"
}

# --- inputs -------------------------------------------------------------
[ -n "${PDF:-}" ] || _cfg_input Supplemental_material.pdf PDF   # Jady et al. supplement
[ -n "${CSV:-}" ] || _cfg_input napRNA_Alu_L1_ACA.csv CSV       # napRNAdb Alu/L1 ACA (hg38)
FASTA="${FASTA:-$OUT/AluACA_HE855917-HE856264.fasta}"           # produced by step 01
# Used by one step each, so they go through _cfg_input like the rest rather than
# assuming the project root -- both actually live in data/ in this checkout.
[ -n "${POLYA_CSV:-}" ] || _cfg_input napRNA_Alu_L1_polyApocketACA.csv POLYA_CSV  # csv_overlap_check.sh
[ -n "${SNODB_TSV:-}" ] || _cfg_input snoDB_All_V2.0.tsv SNODB_TSV                # step 10
FASTA_ID_BASE="${FASTA_ID_BASE:-3000}"                                            # step 08

# --- NCBI ---------------------------------------------------------------
EUTILS="${EUTILS:-https://eutils.ncbi.nlm.nih.gov/entrez/eutils}"
ACC_FROM="${ACC_FROM:-HE855917}"
ACC_TO="${ACC_TO:-HE856264}"
# NCBI allows 3 requests/sec without an API key; steps that loop sleep 0.4s.
NCBI_TOOL="${NCBI_TOOL:-claude_code}"

export PROJ OUT WORK HG38_DIR HG38_FA GENCODE RMSK MAXLEN BIN BEDTOOLS PYTHON \
       PDF CSV POLYA_CSV FASTA SNODB_TSV FASTA_ID_BASE \
       EUTILS ACC_FROM ACC_TO NCBI_TOOL

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

# The project inputs are warnings, not errors: steps 09 and 10 run without
# them. Reported here anyway so a missing PDF surfaces now, with the flag to
# fix it, rather than as a bare "Couldn't open file" from pdftotext in step 02.
for _cfg_pair in "supplemental PDF (steps 02):--pdf:$PDF" \
                 "NapRNAdb CSV (steps 07,08):--csv:$CSV"; do
  _cfg_lbl="${_cfg_pair%%:*}"; _cfg_rest2="${_cfg_pair#*:}"
  _cfg_flag="${_cfg_rest2%%:*}"; _cfg_path="${_cfg_rest2#*:}"
  [ -e "$_cfg_path" ] || echo "WARNING missing $_cfg_lbl: $_cfg_path" \
    "(set it with $_cfg_flag, or put it anywhere under \$PROJ)" >&2
done

[ "$_cfg_missing" -eq 0 ] || {
  echo "" >&2
  echo "Every path is settable by flag or environment variable; see --help." >&2
  exit 1
}
