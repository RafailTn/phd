#!/usr/bin/env bash
# Shared configuration for the chimeric eCLIP pipeline.
# Sourced by every step:  source "$(dirname "$0")/config.sh"
#
# Every path is settable three ways, in increasing precedence:
#   1. built-in default   (derived from the project root, never absolute)
#   2. environment variable
#   3. command-line flag
#
# Because this file is sourced, the flags below work on any step script as well
# as on run_all.sh:
#
#   bash src/chimeric/run_all.sh --species hg19 --source merged SRR30692552
#   bash src/chimeric/run_chimeras.sh --species hg19 --cpus 8 SRR30692552
#   bash src/chimeric/build_indices.sh --dense
#
# Run any script with --help for the full list.
#
# Unlike loci_extraction/config.sh this file does NOT preflight the reference
# files: fetch_refs.sh and build_indices.sh exist precisely to create them, so a
# missing genome here is normal. Each step declares what it needs by calling
# cfg_need itself.
set -euo pipefail

_cfg_usage() {
  cat <<'USAGE'
chimeric eCLIP pipeline - configuration flags

 layout
  --proj DIR         project root (inputs + default output location)
                     default: found by walking up from this file
  --ref DIR          reference directory        [default: $PROJ/ref/chimeric]
  --data DIR         input data directory       [default: $PROJ/data]
  --work DIR         scratch + FASTQ directory  [default: $PROJ/work/chimeric]
  --out DIR          results directory          [default: $PROJ/results/chimeric]

 run
  --species TAG      hg38 | hg19; also the suffix of the genomic outputs
                                                [default: hg38]
  --arm NAME         named arm from the table below; sets --species, --source
                     and the index density together
  --source WHICH     plain | merged | <path>    [default: merged]
                     plain  = snoRNA.txt.fa, the 1936-record catalogue
                     merged = AluACA_snoRNA_merged_nr.fasta, + the 765 AluACAs
  --source-fasta F   the guide catalogue, as an explicit path
  --cpus N           threads                    [default: nproc]
  --srr ACC          sample accession; repeatable. Defaults to the IP and the
                     input control (SRR30692552 SRR30692553).
  --outdir DIR       output directory for a single run_chimeras.sh invocation
                                    [default: $OUT/<arm>/<SRR>]
  --fastq FILE       input reads    [default: $WORK/<SRR>.fastq.gz]

 indices
  --genome-index DIR STAR index for --species   [default: $REF/${SPECIES}_star_index]
  --repeat-index DIR RepBase STAR index         [default: $REF/repbase_star_index]
  --repeat-fa FILE   RepBase consensus FASTA    [default: $REF/repbase/human_repbase.fa]
  --sparse-d N       STAR --genomeSAsparseD. 1 = dense (~32 GB RAM to build,
                     ~29 GB to align); 2 = halved suffix array, ~17 GB.
                                                [default: 1]
  --gen-ram BYTES    STAR --limitGenomeGenerateRAM  [default: 64000000000]
  --dense            shorthand for --sparse-d 1 --gen-ram 64000000000
  --sparse           shorthand for --sparse-d 2 --gen-ram 26000000000
  --sjdb-overhang N  max read length - 1        [default: 139]

 references
  --hg38-dir DIR     directory holding an existing hg38 genome + GENCODE, to
                     reuse instead of downloading  [default: $REF/hg38]
  --genome-fa FILE   primary assembly FASTA (plain or .gz)
  --gencode FILE     annotation GTF (plain or .gz)
  --rmsk FILE        RepeatMasker BED           [default: $REF/rmsk.$SPECIES.bed]
  --guide-bed FILE   guide loci BED             [default: $REF/guide_loci.$SPECIES.bed]
                     Used by run_all.sh for the annotation step. It is NOT
                     passed to sno-chimeras.py unless --source-bed says so.
  --source-bed FILE  hand this BED to sno-chimeras.py as --source_rna_bed, so
                     it drops targets landing back inside the guide's own
                     locus. Off by default: the published runs did not use it,
                     and turning it on changes the chimera counts.

 inputs
  --adapters FILE    second-round adapter FASTA
                                    [default: $REF/se.2.round.adapters.fasta]
  --target-fasta F   target RNA FASTA; repeatable, paired with --target-tag
                                    [default: $DATA/{rRNA,snRNA,tRNA}.fa]
  --target-tag TAG   target RNA tag; repeatable  [default: rRNA snRNA tRNA]
  --alu-fasta FILE   FASTA naming the AluACA records
                                    [default: $DATA/AluACA_union_nr.fasta]
  --published FILE   published hg19 chimeras CSV, for compare_to_published
                              [default: $DATA/DKC1_IP.snoRNA.hg19.chimeras.csv]

 tools
  --bin DIR          directory holding the pipeline's executables
                     default: deps/.pixi/envs/default/bin if present, else PATH
  --python PATH      python interpreter
  --bedtools PATH    bedtools executable

  -h, --help         show this and exit

The reproduction-matrix arms (--arm):

  arm0_hg38_merged          hg38  merged  dense    the headline run
  arm1_hg19_plain           hg19  plain   dense    vs the published 45,810
  arm2_hg38_plain           hg38  plain   dense    catalogue effect
  arm3_hg19_merged          hg19  merged  dense    build-stability of the AluACAs

All arms are dense. --sparse still works for a RAM-limited machine, but names the
run <species>_<source>_sparse so it cannot overwrite a dense result.

Environment variables of the same name in upper snake case (PROJ, REF, DATA,
WORK, OUT, SPECIES, ARM, SOURCE, SOURCE_FASTA, CPUS, SRRS, GENOME_INDEX,
REPEAT_INDEX, REPEAT_FA, SPARSE_D, GEN_RAM, SJDB_OVERHANG, HG38_DIR, GENOME_FA,
GENCODE, RMSK_BED, GUIDE_BED, ADAPTERS, TARGET_FASTA, TARGET_TAGS, ALU_FASTA,
PUBLISHED, BIN, PYTHON, BEDTOOLS) are honoured as defaults; flags win over them.
USAGE
}

# --- flag parsing -------------------------------------------------------
# Unrecognised arguments are left alone so individual steps can define their
# own; only the shared ones are consumed here.
_cfg_rest=()
_cfg_extra=()      # everything after a literal --, passed through to the tool
_cfg_targets=()
_cfg_tags=()
_cfg_srrs=()
while [ "$#" -gt 0 ]; do
  case "$1" in
    --proj)           PROJ="$2";           shift 2 ;;
    --ref)            REF="$2";            shift 2 ;;
    --data)           DATA="$2";           shift 2 ;;
    --work)           WORK="$2";           shift 2 ;;
    --out)            OUT="$2";            shift 2 ;;
    --species)        SPECIES="$2"; _flag_species=1; shift 2 ;;
    --arm)            ARM="$2";            shift 2 ;;
    --source)         SOURCE="$2";  _flag_source=1;  shift 2 ;;
    --source-fasta)   SOURCE_FASTA="$2";   shift 2 ;;
    --cpus)           CPUS="$2";           shift 2 ;;
    --srr)            _cfg_srrs+=("$2");   shift 2 ;;
    --outdir)         OUTDIR="$2";         shift 2 ;;
    --fastq)          FASTQ="$2";          shift 2 ;;
    --source-bed)     SOURCE_BED="$2";     shift 2 ;;
    --genome-index)   GENOME_INDEX="$2";   shift 2 ;;
    --repeat-index)   REPEAT_INDEX="$2";   shift 2 ;;
    --repeat-fa)      REPEAT_FA="$2";      shift 2 ;;
    --sparse-d)       SPARSE_D="$2"; _flag_density=1; shift 2 ;;
    --gen-ram)        GEN_RAM="$2";        shift 2 ;;
    --dense)          SPARSE_D=1; GEN_RAM=${GEN_RAM:-64000000000}; _flag_density=1; shift ;;
    --sparse)         SPARSE_D=2; GEN_RAM=${GEN_RAM:-26000000000}; _flag_density=1; shift ;;
    --sjdb-overhang)  SJDB_OVERHANG="$2";  shift 2 ;;
    --hg38-dir)       HG38_DIR="$2";       shift 2 ;;
    --genome-fa)      GENOME_FA="$2";      shift 2 ;;
    --gencode)        GENCODE="$2";        shift 2 ;;
    --rmsk)           RMSK_BED="$2";       shift 2 ;;
    --guide-bed)      GUIDE_BED="$2";      shift 2 ;;
    --adapters)       ADAPTERS="$2";       shift 2 ;;
    --target-fasta)   _cfg_targets+=("$2"); shift 2 ;;
    --target-tag)     _cfg_tags+=("$2");   shift 2 ;;
    --alu-fasta)      ALU_FASTA="$2";      shift 2 ;;
    --published)      PUBLISHED="$2";      shift 2 ;;
    --bin)            BIN="$2";            shift 2 ;;
    --python)         PYTHON="$2";         shift 2 ;;
    --bedtools)       BEDTOOLS="$2";       shift 2 ;;
    -h|--help)        _cfg_usage; exit 0 ;;
    --)               shift; _cfg_extra+=("$@"); break ;;
    *)                _cfg_rest+=("$1");   shift ;;
  esac
done
set -- ${_cfg_rest[@]+"${_cfg_rest[@]}"}

# --- project layout -----------------------------------------------------
# Walk up from this file until a directory looks like the project root. The
# markers mirror analysis/paths.py: a checkout may be missing the raw inputs but
# still hold the pipeline's own outputs, or vice versa, so anchoring on a single
# file leaves $PROJ pointing at src/.
_cfg_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -z "${PROJ:-}" ]; then
  PROJ="$_cfg_dir"
  while [ "$PROJ" != "/" ]; do
    if [ -d "$PROJ/deps" ] || [ -d "$PROJ/data" ] || [ -e "$PROJ/snoRNA.txt.fa" ] \
       || [ -e "$PROJ/data/snoRNA.txt.fa" ]; then break; fi
    PROJ="$(dirname "$PROJ")"
  done
  # nothing recognisable found: fall back to the grandparent of src/chimeric/
  [ "$PROJ" = "/" ] && PROJ="$(cd "$_cfg_dir/../.." && pwd)"
fi
REF="${REF:-$PROJ/ref/chimeric}"
DATA="${DATA:-$PROJ/data}"
WORK="${WORK:-$PROJ/work/chimeric}"
OUT="${OUT:-$PROJ/results/chimeric}"
SRC="$_cfg_dir"

# --- arm table ----------------------------------------------------------
# One definition of the reproduction matrix, so the driver, build_indices.sh and
# the README cannot drift apart.
#
#            | plain snoRNA.txt.fa   | merged AluACA+snoRNA
#   ---------+-----------------------+----------------------
#   hg38     | arm2                  | arm0 / arm0d
#   hg19     | arm1  <- reproduction | arm3  <- build-stability
# Every arm is dense: that is what produced the committed results, and the
# sparse-vs-dense comparison is settled (see README -- sparse agreed with the
# published run far less well and added calls that were mostly not signal).
# arm0d was "arm0 but dense" and is now the same thing as arm0; it is kept as an
# alias so older commands still land in the right place.
cfg_arm_spec() {   # cfg_arm_spec <arm> -> "<species> <source> <density>"
  case "$1" in
    arm0|arm0_hg38_merged|arm0d|arm0d_hg38_merged_dense)
                                        echo "hg38 merged dense"  ;;
    arm1|arm1_hg19_plain)               echo "hg19 plain  dense"  ;;
    arm2|arm2_hg38_plain)               echo "hg38 plain  dense"  ;;
    arm3|arm3_hg19_merged)              echo "hg19 merged dense"  ;;
    *) return 1 ;;
  esac
}

cfg_arm_name() {   # cfg_arm_name <arm> -> the canonical long name
  case "$1" in
    arm0|arm0_hg38_merged|arm0d|arm0d_hg38_merged_dense)
                                    echo arm0_hg38_merged ;;
    arm1|arm1_hg19_plain)           echo arm1_hg19_plain ;;
    arm2|arm2_hg38_plain)           echo arm2_hg38_plain ;;
    arm3|arm3_hg19_merged)          echo arm3_hg19_merged ;;
    *) return 1 ;;
  esac
}

# --arm sets species, catalogue and density in one go; anything given
# explicitly still wins, since the assignments below are all ${VAR:-...}.
if [ -n "${ARM:-}" ]; then
  _cfg_spec="$(cfg_arm_spec "$ARM")" || {
    echo "unknown arm: $ARM" >&2
    echo "expected one of: arm0 arm0d arm1 arm2 arm3 (or their long names)" >&2
    exit 1; }
  read -r _arm_species _arm_source _arm_density <<<"$_cfg_spec"
  # --arm is a flag, so it outranks an inherited environment variable and only
  # yields to a more specific flag on the same command line. Without this,
  # sourcing config.sh in a shell exports SOURCE/SPECIES and every later
  # --arm silently keeps the old values -- the arm name and the catalogue then
  # disagree, which is invisible until the results are wrong.
  [ -n "${_flag_species:-}" ] || SPECIES="$_arm_species"
  [ -n "${_flag_source:-}" ]  || SOURCE="$_arm_source"
  if [ -z "${_flag_density:-}" ]; then
    if [ "$_arm_density" = dense ]; then
      SPARSE_D=1; GEN_RAM="${GEN_RAM:-64000000000}"
    else
      SPARSE_D=2; GEN_RAM="${GEN_RAM:-26000000000}"
    fi
  fi
  ARM="$(cfg_arm_name "$ARM")"
fi

# The reverse lookup: results always live under $OUT/<arm>/<SRR>, so a run
# configured by --species/--source/--dense still lands in the arm directory it
# corresponds to. A combination outside the table gets a descriptive name.
cfg_arm_for() {   # cfg_arm_for <species> <source> <density> -> arm name
  case "$1 $2 $3" in
    "hg38 merged dense")  echo arm0_hg38_merged ;;
    "hg19 plain  dense"|"hg19 plain dense")   echo arm1_hg19_plain ;;
    "hg38 plain  dense"|"hg38 plain dense")   echo arm2_hg38_plain ;;
    "hg19 merged dense")  echo arm3_hg19_merged ;;
    *) echo "${1}_${2}_${3}" ;;
  esac
}

# --- run parameters -----------------------------------------------------
SPECIES="${SPECIES:-hg38}"
CPUS="${CPUS:-$(nproc)}"
# Dense by default: it is what produced every committed result, and the named
# arms all specify it. A RAM-limited machine passes --sparse, which names the run
# so it cannot be mistaken for, or overwrite, a dense one.
SPARSE_D="${SPARSE_D:-1}"
GEN_RAM="${GEN_RAM:-64000000000}"
# sjdbOverhang = max read length - 1. Reads are 150 nt minus a 10 nt UMI, so 139.
SJDB_OVERHANG="${SJDB_OVERHANG:-139}"

# The IP and its input control. See data/ for the accession-to-GSM mapping.
if [ "${#_cfg_srrs[@]}" -gt 0 ]; then SRRS="${_cfg_srrs[*]}"; fi
SRRS="${SRRS:-SRR30692552 SRR30692553}"
IP="${IP:-${SRRS%% *}}"

# --- input lookup -------------------------------------------------------
# An input is looked for in the directories it conventionally lives in, then
# anywhere under $PROJ. The search never leaves the project: the current
# directory is deliberately not consulted, so running from somewhere else
# cannot silently pick up a same-named file. Flags override the lookup entirely.
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
  for d in "$@" "$DATA" "$PROJ"; do
    if [ -e "$d/$base" ]; then printf -v "$var" '%s' "$d/$base"; return 0; fi
  done
  _cfg_build_index
  # Shallowest match wins, ties broken alphabetically, so which file is chosen
  # never depends on the order the filesystem happened to return.
  hit=$(printf '%s\n' "$_cfg_index" | awk -v b="$base" '
          { n = split($0, p, "/"); if (p[n] == b) print n "\t" $0 }' \
        | sort -k1,1n -k2,2 | head -1 | cut -f2-)
  # Report the conventional path when missing, so the error names where it goes.
  printf -v "$var" '%s' "${hit:-$DATA/$base}"
}

# --- inputs -------------------------------------------------------------
# The guide ("source") catalogue. Keeping the real snoRNAs in alongside the
# AluACAs is deliberate: bowtie2 scores a read against the whole catalogue at
# once, so an AluACA call means the AluACA beat every snoRNA, rather than being
# the only thing on offer. AluACA records are the ones with a .id3xxx suffix.
SOURCE="${SOURCE:-merged}"
if [ -z "${SOURCE_FASTA:-}" ]; then
  case "$SOURCE" in
    plain)  _cfg_input snoRNA.txt.fa SOURCE_FASTA ;;
    merged) _cfg_input AluACA_snoRNA_merged_nr.fasta SOURCE_FASTA ;;
    *)      SOURCE_FASTA="$SOURCE" ;;   # an explicit path
  esac
fi

# These two normally live under ref/ rather than data/, so it is tried first.
[ -n "${ADAPTERS:-}" ]  || _cfg_input se.2.round.adapters.fasta ADAPTERS "$REF"
[ -n "${ALU_FASTA:-}" ] || _cfg_input AluACA_union_nr.fasta ALU_FASTA
[ -n "${PUBLISHED:-}" ] || _cfg_input DKC1_IP.snoRNA.hg19.chimeras.csv PUBLISHED

# Target catalogues the guide arm is mapped against, tags in the same order.
if [ "${#_cfg_targets[@]}" -gt 0 ]; then TARGET_FASTA="${_cfg_targets[*]}"; fi
if [ "${#_cfg_tags[@]}" -gt 0 ];    then TARGET_TAGS="${_cfg_tags[*]}"; fi
if [ -z "${TARGET_FASTA:-}" ]; then
  _cfg_input rRNA.fa  _cfg_t_r
  _cfg_input snRNA.fa _cfg_t_sn
  _cfg_input tRNA.fa  _cfg_t_t
  TARGET_FASTA="$_cfg_t_r $_cfg_t_sn $_cfg_t_t"
fi
TARGET_TAGS="${TARGET_TAGS:-rRNA snRNA tRNA}"
# The tag treated as rRNA when deciding which chimeras are ribosomal.
TARGET_RRNA_TAG="${TARGET_RRNA_TAG:-rRNA}"

# --- reference data -----------------------------------------------------
# One statement of which files each build is made of. fetch_refs.sh downloads
# exactly these names, build_indices.sh builds from them and run_all.sh
# annotates with them, so the three cannot drift apart.
GENCODE_RELEASE="${GENCODE_RELEASE:-47}"

# The only statement of which files a build is made of. fetch_refs.sh downloads
# these names, build_indices.sh builds from them and run_all.sh annotates with
# them, so the three cannot drift apart.
#
# hg19 uses the vNNlift37 mapping, so the annotation release is the same across
# builds and only the coordinates differ.
cfg_ref_names() {   # cfg_ref_names <build> -> "<fasta> <gtf> <gencode subdir|->"
  case "$1" in
    hg38) echo "GRCh38.primary_assembly.genome.fa" \
               "gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz" "-" ;;
    hg19) echo "GRCh37.primary_assembly.genome.fa" \
               "gencode.v${GENCODE_RELEASE}lift37.annotation.gtf.gz" "GRCh37_mapping/" ;;
    *) return 1 ;;
  esac
}

if _cfg_names="$(cfg_ref_names "$SPECIES")"; then
  read -r GENOME_FA_NAME GENCODE_NAME GENCODE_SUBDIR <<<"$_cfg_names"
  [ "$GENCODE_SUBDIR" = "-" ] && GENCODE_SUBDIR=
else
  # An unknown build is allowed, but then --genome-fa and --gencode are required.
  GENOME_FA_NAME= ; GENCODE_NAME= ; GENCODE_SUBDIR=
fi

# An existing genome outside the repo is reused rather than re-downloaded; the
# files are the same GENCODE release either way. $HOME/Downloads/hg38 matches
# the default loci_extraction/config.sh already uses.
HG38_DIR="${HG38_DIR:-$REF/hg38}"
_cfg_ref_dir="$REF/$SPECIES"
if [ "$SPECIES" = hg38 ]; then
  for _d in "$HG38_DIR" "$REF/hg38" "$HOME/Downloads/hg38"; do
    if [ -f "$_d/$GENOME_FA_NAME" ]; then _cfg_ref_dir="$_d"; break; fi
  done
fi
REF_SPECIES_DIR="${REF_SPECIES_DIR:-$_cfg_ref_dir}"

if [ -z "${GENOME_FA:-}" ] && [ -n "$GENOME_FA_NAME" ]; then
  GENOME_FA="$REF_SPECIES_DIR/$GENOME_FA_NAME"
fi
if [ -z "${GENCODE:-}" ] && [ -n "$GENCODE_NAME" ]; then
  GENCODE="$REF_SPECIES_DIR/$GENCODE_NAME"
fi
GENOME_FA="${GENOME_FA:-}"
GENCODE="${GENCODE:-}"
RMSK_BED="${RMSK_BED:-$REF/rmsk.$SPECIES.bed}"
GUIDE_BED="${GUIDE_BED:-$REF/guide_loci.$SPECIES.bed}"

GENOME_INDEX="${GENOME_INDEX:-$REF/${SPECIES}_star_index}"
REPEAT_INDEX="${REPEAT_INDEX:-$REF/repbase_star_index}"
[ -n "${REPEAT_FA:-}" ] || _cfg_input human_repbase.fa REPEAT_FA "$REF/repbase" "$REF"

# Download endpoints, for fetch_refs.sh.
GENCODE_URL="${GENCODE_URL:-https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}}"
UCSC_URL="${UCSC_URL:-https://hgdownload.soe.ucsc.edu/goldenPath}"

# Name the arm from the resolved configuration when it was not given directly,
# so every run has one and results always land under $OUT/<arm>/<SRR>.
if [ -z "${ARM:-}" ]; then
  if [ "$SPARSE_D" = 1 ]; then _cfg_density=dense; else _cfg_density=sparse; fi
  ARM="$(cfg_arm_for "$SPECIES" "$SOURCE" "$_cfg_density")"
fi

# --- tools --------------------------------------------------------------
# Prefer the project's pixi env when it exists; otherwise fall back to PATH so
# the pipeline runs under a module system, conda env or plain container.
BIN="${BIN:-$PROJ/deps/.pixi/envs/default/bin}"
[ -d "$BIN" ] && case ":$PATH:" in *":$BIN:"*) ;; *) PATH="$BIN:$PATH" ;; esac
# The bin/ wrappers are how sno-chimeras.py reaches its helper scripts: it
# invokes them as bare command names.
case ":$PATH:" in *":$SRC/bin:"*) ;; *) PATH="$SRC/bin:$PATH" ;; esac
export PATH

if [ -z "${PYTHON:-}" ]; then
  if [ -x "$BIN/python" ];    then PYTHON="$BIN/python"
  elif [ -x "$BIN/python3" ]; then PYTHON="$BIN/python3"
  else PYTHON="$(command -v python3 || true)"; fi
fi
if [ -z "${BEDTOOLS:-}" ]; then
  if [ -x "$BIN/bedtools" ]; then BEDTOOLS="$BIN/bedtools"
  else BEDTOOLS="$(command -v bedtools || true)"; fi
fi

# sno-chimeras.py resolves its built-in adapter FASTAs relative to this.
CHIMERIC_REF_DIR="${CHIMERIC_REF_DIR:-$REF}"

# Single-run overrides, empty unless given. OUTDIR and FASTQ default per sample
# in run_chimeras.sh; SOURCE_BED stays empty so the pipeline behaves as it did.
OUTDIR="${OUTDIR:-}"
FASTQ="${FASTQ:-}"
SOURCE_BED="${SOURCE_BED:-}"

export PROJ SRC REF DATA WORK OUT SPECIES ARM SOURCE SOURCE_FASTA CPUS SRRS IP \
       GENOME_INDEX REPEAT_INDEX REPEAT_FA SPARSE_D GEN_RAM SJDB_OVERHANG \
       HG38_DIR REF_SPECIES_DIR GENOME_FA GENCODE RMSK_BED GUIDE_BED \
       GENOME_FA_NAME GENCODE_NAME GENCODE_SUBDIR GENCODE_RELEASE \
       GENCODE_URL UCSC_URL \
       ADAPTERS TARGET_FASTA TARGET_TAGS TARGET_RRNA_TAG ALU_FASTA PUBLISHED \
       BIN PYTHON BEDTOOLS CHIMERIC_REF_DIR OUTDIR FASTQ SOURCE_BED

# --- preflight helpers --------------------------------------------------
# Offered rather than run: fetch_refs.sh and build_indices.sh exist to create
# the very files a blanket check would demand, so each step says what it needs.
_cfg_missing=0
cfg_need() {   # cfg_need <label> <path> <how-to-set>
  if [ -z "${2:-}" ]; then
    echo "MISSING $1: could not be located automatically  ($3)" >&2; _cfg_missing=1
  elif [ ! -e "$2" ]; then
    echo "MISSING $1: $2  ($3)" >&2; _cfg_missing=1
  fi
}
cfg_check() {  # cfg_check -- exit non-zero if any cfg_need failed
  [ "$_cfg_missing" -eq 0 ] && return 0
  echo "" >&2
  echo "Every path is settable by flag or environment variable; see --help." >&2
  exit 1
}
