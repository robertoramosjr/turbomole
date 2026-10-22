#!/usr/bin/env python3
"""Prepare xtb (GFN2, --hess) G_thermal runs for the water-ensemble W4 conformers.

Every conformer with a W4 single point (<root>/w4_sp/<state>/<cNNNN>/W4_DONE) gets its TURBOMOLE
`coord` (bohr) converted to <root>/xtb_gtherm_w4/<state>/<cNNNN>/{mol.xyz,charge.txt}; the absolute
dirs are written to a joblist for job_xtb_hess.sh.  The Hessian is taken at the same geometry as the
W4 energy (gas-phase xtb, as in the validated pilot).

Usage:
  python prepare_xtb_gtherm_w4.py --root ~/work_turbomole/artepillin_C_studies --joblist xtb_gtherm_w4/joblist.txt
"""
import argparse
import os

from prepare_xtb_gtherm_pilot import STATES, coord_to_xyz


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--joblist", required=True, help="relative to --root; overwritten")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    dirs = []
    for state, chg in STATES.items():
        base = os.path.join(root, "w4_sp", state)
        for c in sorted(os.listdir(base)):
            src = os.path.join(base, c, "coord")
            if not os.path.isfile(os.path.join(base, c, "W4_DONE")) or not os.path.isfile(src):
                print(f"skip {state}/{c}: no W4_DONE/coord")
                continue
            dst = os.path.join(root, "xtb_gtherm_w4", state, c)
            os.makedirs(dst, exist_ok=True)
            with open(os.path.join(dst, "mol.xyz"), "w") as out:
                out.write(coord_to_xyz(src))
            with open(os.path.join(dst, "charge.txt"), "w") as out:
                out.write(f"{chg}\n")
            dirs.append(dst)
    os.makedirs(os.path.dirname(os.path.join(root, a.joblist)), exist_ok=True)
    with open(os.path.join(root, a.joblist), "w") as fh:
        fh.write("".join(d + "\n" for d in dirs))
    print(f"prepared {len(dirs)} dirs")


if __name__ == "__main__":
    main()
