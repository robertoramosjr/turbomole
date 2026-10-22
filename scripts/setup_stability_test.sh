#!/bin/bash
# Prepare wavefunction-stability test dirs (SKILL turbomole-pipeline, hard rule 7; manual 8.4.5).
# Copies only the INPUTS of a converged calc (no twoint/gradient/hessian) into
# <outroot>/<tag>/ and rewrites $scfinstab/$soes. Usage:
#   setup_stability_test.sh <converged_dir> <tag> <singlet|triplet|urpa> [outroot]
set -eu
SRC=$(readlink -f "$1"); TAG=$2; KIND=$3; OUT=${4:-$HOME/work_turbomole/metal_complexes/stability}
[ -f "$SRC/GEO_OPT_CONVERGED" ] || { echo "no GEO_OPT_CONVERGED in $SRC"; exit 1; }
D=$OUT/$TAG; mkdir -p "$D"
for f in control coord basis auxbasis mos alpha beta; do [ -f "$SRC/$f" ] && cp "$SRC/$f" "$D/"; done
cd "$D"
# drop old excited-state / instability blocks, keep everything else (functional, grid, $disp3, $scfconv, ...)
python3 - "$KIND" <<'PY'
import re, sys
kind = sys.argv[1]
lines = open("control").read().split("\n")
out, skip = [], False
for l in lines:
    if l.startswith("$"):
        skip = l.split()[0] in ("$scfinstab", "$soes", "$start")
        if l.startswith("$end"):
            continue
    if not skip:
        out.append(l)
new = [f"$scfinstab {kind}", "$soes", " all 1"]
if kind in ("singlet", "triplet"):
    new.append("$start vector generation")  # BS start MOs, only written if a totally-symmetric triplet instability exists
open("control", "w").write("\n".join(out).rstrip("\n") + "\n" + "\n".join(new) + "\n$end\n")
PY
cp "$SRC/control" control.orig
echo "prepared $D ($KIND)"
