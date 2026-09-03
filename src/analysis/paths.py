"""Input/output/tool lookup shared by the analysis scripts.

Mirrors loci_extraction/config.sh so the two halves of the project agree on
where things live:

  inputs   ./  ->  $PROJ/  ->  $PROJ/data/          (data/ is the server layout)
  outputs  $OUT, or $PROJ when OUT is unset
  tools    $BIN -> the pixi env under $PROJ/deps -> PATH

$PROJ is taken from the environment when set, otherwise found by walking up
from this file until a directory looks like the project root.
"""

import os
import shutil

_HERE = os.path.dirname(os.path.abspath(__file__))

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
    """The directories find_input consults, in order."""
    return [os.curdir, proj(), os.path.join(proj(), "data")]


def find_input(basename):
    """First search directory holding the file; else the $PROJ path, so the
    caller can report a sensible location when it is missing."""
    for d in search_dirs():
        cand = os.path.join(d, basename)
        if os.path.exists(cand):
            return os.path.normpath(cand)
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
    one, naming the directories that were searched."""
    import sys
    for label, path in pairs:
        if not path or not os.path.exists(path):
            sys.exit(f"missing {label}: {path}\n  looked in "
                     + ", ".join(os.path.abspath(d) for d in search_dirs()))
