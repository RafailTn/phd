#!/usr/bin/env bash
# Pack everything another machine needs, except the FASTQs and the STAR indices
# (fetch_fastq.sh downloads the first, build_indices.sh builds the second --
# that is the point of the exercise).
#
#   bash src/chimeric/make_bundle.sh [outfile.tar.gz]
#
# On the other machine:
#   mkdir phd && tar -xzf chimeric-bundle.tar.gz -C phd && cd phd
#   pixi install --manifest-path deps/pixi.toml
#   bash src/chimeric/run_all.sh arm1 arm3
#
# No `export PROJ=$PWD` is needed: config.sh finds the project root by walking
# up from itself, and the bin/ wrappers are self-locating.
#
# The file list is derived from config.sh rather than restated here, so an input
# added to the pipeline cannot be forgotten by the bundle.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

OUT_TAR=${1:-$WORK/chimeric-bundle.tar.gz}

# Both guide catalogues, whichever arm the other machine ends up running.
PLAIN=$(SOURCE=plain  SOURCE_FASTA= bash -c 'source "$0" >/dev/null; echo "$SOURCE_FASTA"' "$SRC/config.sh")
MERGED=$(SOURCE=merged SOURCE_FASTA= bash -c 'source "$0" >/dev/null; echo "$SOURCE_FASTA"' "$SRC/config.sh")

FILES=(
  src/chimeric
  src/paths.py                # the chimeric scripts import it from src/
  deps/pixi.toml deps/pixi.lock
)
# Everything config.sh points at, as a path relative to the project root.
for f in "$PLAIN" "$MERGED" "$ALU_FASTA" $TARGET_FASTA "$PUBLISHED" \
         "$ADAPTERS" "$REPEAT_FA" "$REF/guide_loci.hg38.bed"; do
  case "$f" in
    "$PROJ"/*) FILES+=("${f#"$PROJ"/}") ;;
    *) echo "outside the project, not bundled: $f" >&2 ;;
  esac
done

cd "$PROJ"
missing=0
for f in "${FILES[@]}"; do
  [ -e "$f" ] || { echo "missing: $f" >&2; missing=1; }
done
[ "$missing" -eq 0 ] || { echo "" >&2; echo "nothing was written" >&2; exit 1; }

mkdir -p "$(dirname "$OUT_TAR")"
tar -czf "$OUT_TAR" --exclude='__pycache__' "${FILES[@]}"
echo "wrote $OUT_TAR ($(du -h "$OUT_TAR" | cut -f1))"
tar -tzf "$OUT_TAR" | wc -l | xargs echo "entries:"
