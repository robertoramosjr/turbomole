#!/usr/bin/env python3
"""Prepare the xtb G_thermal validation pilot (GFN2-xTB --bhess vs HSE06 freeh).

For every conformer with a HSE06 `freeh` correction (layer_b/freeh_thermal_correction.csv) copy the
HSE06-optimised geometry (layer_b/labeled_jobs/<job>/coord, TURBOMOLE bohr) as an xyz (angstrom) into
<root>/xtb_gtherm_pilot/<state>/<job>/mol.xyz and append the dir to a joblist.

Usage:
  python prepare_xtb_gtherm_pilot.py --root ~/work_turbomole/artepillin_C_studies --joblist xtb_gtherm_pilot/joblist.txt
"""
import argparse
import csv
import os

BOHR_TO_ANG = 0.529177210903
STATES = {"neutral": 0, "monoanion_a": -1, "monoanion_b": -1, "dianion": -2}


def coord_to_xyz(path):
    atoms = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("$coord"):
                continue
            if ln.startswith("$"):
                break
            p = ln.split()
            if len(p) >= 4:
                atoms.append((p[3].capitalize(), [float(v) * BOHR_TO_ANG for v in p[:3]]))
    lines = [str(len(atoms)), "from TURBOMOLE coord (bohr->angstrom)"]
    lines += [f"{s:2s} {x:16.8f} {y:16.8f} {z:16.8f}" for s, (x, y, z) in atoms]
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--joblist", required=True, help="relative to --root; overwritten")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    dirs = []
    for state, chg in STATES.items():
        base = os.path.join(root, f"artepillin_c_{state}_d3_hse06", "layer_b")
        with open(os.path.join(base, "freeh_thermal_correction.csv")) as fh:
            jobs = [r["job_dir"] for r in csv.DictReader(fh)]
        for j in jobs:
            src = os.path.join(base, "labeled_jobs", j, "coord")
            if not os.path.isfile(src):
                print(f"skip {state}/{j}: no coord")
                continue
            dst = os.path.join(root, "xtb_gtherm_pilot", state, j)
            os.makedirs(dst, exist_ok=True)
            with open(os.path.join(dst, "mol.xyz"), "w") as out:
                out.write(coord_to_xyz(src))
            with open(os.path.join(dst, "charge.txt"), "w") as out:
                out.write(f"{chg}\n")
            dirs.append(dst)
    with open(os.path.join(root, a.joblist), "w") as fh:
        fh.write("".join(d + "\n" for d in dirs))
    print(f"prepared {len(dirs)} dirs")


if __name__ == "__main__":
    main()
