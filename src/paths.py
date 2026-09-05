"""Input/output/tool lookup shared by the analysis and chimeric scripts.

Mirrors loci_extraction/config.sh and chimeric/config.sh so the halves of the
project agree on where things live:

  inputs   $PROJ/  ->  $PROJ/data/  ->  anywhere under $PROJ
  outputs  $OUT, or $PROJ when OUT is unset
  tools    $BIN -> the pixi env under $PROJ/deps -> PATH

$PROJ is taken from the environment when set, otherwise found by walking up
from this file until a directory looks like the project root.

Input lookup never leaves $PROJ. The current directory is deliberately not
searched: it is usually the project root anyway, and when it is not, a
same-named file in an unrelated directory would be picked up silently.
"""

import os
import shutil

_HERE = os.path.dirname(os.path.abspath(__file__))

# Skipped by the recursive input search. deps/ is a vendored conda environment
# of ~28k files that are not ours to match against; the STAR indices are large
# and hold no inputs. Both would only slow the walk down and invite collisions.
_PRUNE_DIRS = {".git", "deps", "__pycache__"}
_PRUNE_SUFFIXES = ("_star_index",)

# Any one of these at a directory's root -- or inside its data/ -- marks it as
# the project root.  Several markers rather than one because a checkout may be
# missing the raw inputs but still hold the pipeline's own outputs, or vice
# versa; anchoring on a single file leaves $PROJ pointing at scripts/.
_MARKERS = (
    "napRNA_Alu_L1_ACA.csv",
    "snoRNA.txt.fa",
    "AluACA_union_nr.fasta",
    "Supplemental_material.pdf",
)


def proj():
    """The project root."""
    env = os.environ.get("PROJ")
    if env:
        return env
    p = _HERE
    while p != os.path.dirname(p):
        if os.path.isdir(os.path.join(p, "deps")) \
           or any(os.path.exists(os.path.join(p, m)) for m in _MARKERS) \
           or any(os.path.exists(os.path.join(p, "data", m)) for m in _MARKERS):
            return p
        p = os.path.dirname(p)
    return os.path.dirname(_HERE)


def search_dirs():
    """The directories find_input checks before falling back to a full walk."""
    return [proj(), os.path.join(proj(), "data")]


_index = None


def _repo_index():
    """basename -> [paths], for every file under $PROJ outside the pruned dirs.

    Built once on the first lookup that misses the conventional directories, so
    a checkout with everything in data/ never pays for the walk.
    """
    global _index
    if _index is None:
        _index = {}
        # followlinks stays False: a symlink out of the tree must not become a
        # way for the search to leave $PROJ.
        for dirpath, dirnames, filenames in os.walk(proj()):
            dirnames[:] = [d for d in dirnames
                           if d not in _PRUNE_DIRS
                           and not d.endswith(_PRUNE_SUFFIXES)]
            for fn in filenames:
                _index.setdefault(fn, []).append(os.path.join(dirpath, fn))
    return _index


def find_input(basename):
    """Locate an input by name, without leaving $PROJ.

    The conventional homes are checked first, then anywhere under the project
    root. Of several matches the shallowest wins, ties broken alphabetically,
    so the answer does not depend on filesystem order. Returns the $PROJ path
    when nothing is found, so the caller can report a sensible location.
    """
    for d in search_dirs():
        cand = os.path.join(d, basename)
        if os.path.exists(cand):
            return os.path.normpath(cand)
    hits = _repo_index().get(basename)
    if hits:
        return min(hits, key=lambda p: (p.count(os.sep), p))
    return os.path.join(proj(), basename)


def out_path(basename):
    """Deliverables go to $OUT, defaulting to the project root."""
    return os.path.join(os.environ.get("OUT") or proj(), basename)


def find_tool(name, explicit=None):
    """Locate an executable: explicit path, then $NAME, then $BIN, then the
    project's pixi env, then PATH.  Returns None if nothing is found."""
    if explicit:
        return explicit
    env = os.environ.get(name.upper())
    if env:
        return env
    for d in (os.environ.get("BIN"),
              os.path.join(proj(), "deps", ".pixi", "envs", "default", "bin")):
        if d:
            cand = os.path.join(d, name)
            if os.access(cand, os.X_OK):
                return cand
    return shutil.which(name)


def require(pairs):
    """pairs: (label, path).  Exit with a useful message on the first missing
    one, naming where it was looked for."""
    import sys
    for label, path in pairs:
        if not path or not os.path.exists(path):
            sys.exit(f"missing {label}: {path}\n  looked in "
                     + ", ".join(os.path.abspath(d) for d in search_dirs())
                     + f", and anywhere under {proj()}")
