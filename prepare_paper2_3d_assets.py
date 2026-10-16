#!/usr/bin/env python3
"""
prepare_paper2_3d_assets.py

Collects the 3D inputs for the Paper 2 renders (VESTA/Jmol): for each
protonation microstate, the r2SCAN-3c geometry of the most populated conformer
(single .xyz), a multi-frame .xyz of the top-N conformers by r2SCAN-3c
population (frame comment = job, rank, population), and any existing Gaussian
cube files of the Paper 1 reference conformer (transition density), copied unchanged.
NTO/transition-density cubes for the Paper 2 microstates do NOT exist yet (would need new calculations).

Usage:
    python prepare_paper2_3d_assets.py --project-root ~/work_turbomole/artepillin_C_studies \
        --outdir ~/work_turbomole/artepillin_C_studies/paper_2/Figures/3D_for_VESTA --top 10
"""
import argparse
import csv
import os
import shutil

BOHR = 0.529177210903
STATES = {"neutral": "N", "monoanion_a": "A", "monoanion_b": "B", "dianion": "D"}


def read_coord(path):
    atoms, on = [], False
    with open(path) as f:
        for ln in f:
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


def write_xyz(fh, atoms, comment):
    fh.write(f"{len(atoms)}\n{comment}\n")
    for el, x, y, z in atoms:
        fh.write(f"{el:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")


def main():
    ap = argparse.ArgumentParser(description="Collect 3D assets for Paper 2 renders.")
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--top", type=int, default=10)
    a = ap.parse_args()
    root, out = os.path.expanduser(a.project_root), os.path.expanduser(a.outdir)
    os.makedirs(out, exist_ok=True)
    for state, tag in STATES.items():
        d = os.path.join(root, f"artepillin_c_{state}_d3_hse06", "layer_a")
        with open(os.path.join(d, "population_table_dft.csv")) as f:
            rows = sorted(csv.DictReader(f), key=lambda r: -float(r["r2scan3c_population"]))
        pop_tot = sum(float(r["r2scan3c_population"]) for r in rows)
        frames = []
        for rank, r in enumerate(rows[: a.top], 1):
            atoms = read_coord(os.path.join(d, "conformer_jobs", r["job_dir"], "coord"))
            frames.append((rank, r, atoms))
        r0, at0 = frames[0][1], frames[0][2]
        with open(os.path.join(out, f"{tag}_{state}_top1.xyz"), "w") as fh:
            write_xyz(fh, at0, f"{state} {r0['job_dir']} r2SCAN-3c pop {float(r0['r2scan3c_population']) / pop_tot:.4f}")
        with open(os.path.join(out, f"{tag}_{state}_top{a.top}_ensemble.xyz"), "w") as fh:
            for rank, r, atoms in frames:
                write_xyz(fh, atoms, f"{state} rank {rank} {r['job_dir']} pop {float(r['r2scan3c_population']) / pop_tot:.4f} "
                                      f"rel_E_r2SCAN3c {float(r['r2scan3c_relative_energy_kcal']):.3f} kcal/mol")
        print(state, "top1 =", r0["job_dir"], "natoms", len(at0))
    cubdir = os.path.join(out, "cube_files")
    os.makedirs(cubdir, exist_ok=True)
    ref = os.path.join(root, "artepillin_C_d3_hse06", "td.cub")
    if os.path.isfile(ref):
        shutil.copy(ref, os.path.join(cubdir, "N_paper1_reference_conformer_td.cub"))
    print("cubes copied:", len(os.listdir(cubdir)))


if __name__ == "__main__":
    main()
