#!/usr/bin/env python3
"""
water_reopt_conformer.py

Build ONE directory for r2SCAN-3c/def2-mTZVP + COSMO(water) geometry re-optimisation of a CREST/GFN2-ALPB conformer:
copies control/basis/auxbasis/mos from a COSMO-ready template (see water_tier2_sp.py template), strips stale
optimisation bookkeeping, sets $scfiterlimit and writes coord from crest_conformers.xyz (1-based index).

Usage:
  python water_reopt_conformer.py --xyz water_ensemble/dianion/crest_conformers.xyz --template tier2/dianion/template \
         --index 77 --outdir reopt/dianion/c0077
"""
import argparse
import os
import re
import shutil

BOHR = 1.0 / 0.529177210903


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--template", required=True)
    ap.add_argument("--index", type=int, required=True, help="1-based conformer index in the xyz")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--scfiterlimit", type=int, default=150)
    a = ap.parse_args()

    lines = open(a.xyz).read().split("\n")
    nat = int(lines[0])
    i0 = (a.index - 1) * (nat + 2)
    atoms = [l.split() for l in lines[i0 + 2:i0 + 2 + nat]]
    if len(atoms) != nat or any(len(x) < 4 for x in atoms):
        raise SystemExit(f"conformer {a.index} not found in {a.xyz}")
    os.makedirs(a.outdir, exist_ok=True)
    for f in ("control", "basis", "auxbasis", "mos"):
        shutil.copy(os.path.join(a.template, f), a.outdir)
    p = os.path.join(a.outdir, "control")
    c = open(p).read()
    c = re.sub(r"^\$(scfiterlimit|last step|actual step).*\n", "", c, flags=re.M)
    c = c.replace("$end", f"$scfiterlimit {a.scfiterlimit}\n$end", 1)
    open(p, "w").write(c)
    with open(os.path.join(a.outdir, "coord"), "w") as f:
        f.write("$coord\n")
        for at in atoms:
            f.write(f"{float(at[1]) * BOHR:20.12f}{float(at[2]) * BOHR:20.12f}{float(at[3]) * BOHR:20.12f}     {at[0].lower()}\n")
        f.write("$end\n")


if __name__ == "__main__":
    main()
