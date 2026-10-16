#!/usr/bin/env python3
"""
xtb_solvation_on_dft_geometries.py

Apples-to-apples check of the conformer-dependent SOLVATION energy: runs GFN2-xTB
single points (gas and ALPB water) on exactly the r2SCAN-3c geometries that the
HSE06/COSMO(OC) numbers were computed on, so that
    dG_solv(ALPB, conformer) = E_ALPB - E_gas        (GFN2)
can be compared conformer by conformer with
    dG_solv(COSMO)  = E_COSMO(OC) - E_gas(HSE06)     (from analyze_vacuum_vs_water_pka.py).
Needed because the earlier GFN2/ALPB reranking re-optimised in water (different
geometry) and could not be compared directly with the DFT/COSMO energies.

Conformer sets (same as analyze_vacuum_vs_water_pka.py): labelled subset of every
state + the water-guided/top-30 extra conformers.

Usage (run through sbatch, ~0.2 core-h total):
    python xtb_solvation_on_dft_geometries.py --project-root ~/work_turbomole/artepillin_C_studies \
        --xtb ~/software/TURBOMOLE/bin/em64t-unknown-linux-gnu/xtb --workdir xtb_solv_check \
        --output xtb_solvation_on_dft_geometries.csv
"""
import argparse
import csv
import os
import re
import subprocess

BOHR = 0.529177210903
STATES = {"N": "neutral", "A": "monoanion_a", "B": "monoanion_b", "D": "dianion"}
EXT_DIR = {"A": "top30_ext", "B": "water_ext", "D": "water_ext"}


def read_coord(path):
    atoms, on = [], False
    for ln in open(path):
        if ln.startswith("$coord"):
            on = True
            continue
        if on:
            if ln.startswith("$"):
                break
            p = ln.split()
            if len(p) >= 4:
                atoms.append((p[3].capitalize(), *(float(v) * BOHR for v in p[:3])))
    return atoms


def write_xyz(path, atoms):
    with open(path, "w") as f:
        f.write(f"{len(atoms)}\n\n")
        for el, x, y, z in atoms:
            f.write(f"{el:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")


def run_xtb(xtb, xyz, chrg, alpb, cwd):
    cmd = [xtb, xyz, "--gfn", "2", "--sp", "--chrg", str(chrg)] + (["--alpb", "water"] if alpb else [])
    out = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True).stdout
    m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", out)
    return float(m[-1]) if m else None


def jobs_for(root, s):
    base = os.path.join(root, f"artepillin_c_{STATES[s]}_d3_hse06")
    res = []
    rows = list(csv.DictReader(open(os.path.join(base, "layer_b", "cosmo_energy.csv"))))
    for r in rows:
        res.append(("old", r["job_dir"], os.path.join(base, "layer_a", "conformer_jobs", r["job_dir"], "coord")))
    ed = os.path.join(base, "layer_b", EXT_DIR.get(s, "none"))
    if s in EXT_DIR and os.path.exists(os.path.join(ed, "cosmo_energy.csv")):
        for r in csv.DictReader(open(os.path.join(ed, "cosmo_energy.csv"))):
            j = r["job_dir"]
            if s == "A":
                res.append(("new", j, os.path.join(base, "layer_a", "conformer_jobs", j, "coord")))
            else:
                res.append(("new", j, os.path.join(base, "layer_a", "water_ext", "conformer_jobs", j, "coord")))
    return base, res


def main():
    ap = argparse.ArgumentParser(description="GFN2 gas/ALPB single points on the DFT geometries.")
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--xtb", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    root = os.path.expanduser(a.project_root)
    os.makedirs(a.workdir, exist_ok=True)
    out = []
    for s in STATES:
        base, jobs = jobs_for(root, s)
        chrg = int(open(os.path.join(base, "CHARGE.txt")).read().split()[0])
        for tag, j, coord in jobs:
            if not os.path.exists(coord):
                print("missing", coord)
                continue
            d = os.path.join(a.workdir, f"{s}_{tag}_{j}")
            os.makedirs(d, exist_ok=True)
            write_xyz(os.path.join(d, "g.xyz"), read_coord(coord))
            eg = run_xtb(a.xtb, "g.xyz", chrg, False, d)
            ew = run_xtb(a.xtb, "g.xyz", chrg, True, d)
            out.append([s, tag, j, eg, ew, (ew - eg) * 627.509474 if eg and ew else ""])
            print(s, tag, j, eg, ew, flush=True)
    with open(a.output, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["state", "set", "job", "E_gas_xtb_Eh", "E_alpb_xtb_Eh", "dG_solv_alpb_kcal"])
        w.writerows(out)


if __name__ == "__main__":
    main()
