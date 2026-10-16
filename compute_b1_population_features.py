#!/usr/bin/env python3
"""
compute_b1_population_features.py

Extracts bond order (Mulliken shared electron number, SEN) for the
specific C=O bond and Bader net charges for the specific carbonyl C
and O atoms -- for every conformer in a job-directory table, from the
population-analysis outputs of job_population_analysis.sh (pop.out for
SEN, ACF.dat for Bader). Carbonyl atom identification reuses the same
geometry-based detection as identify_co_stretch_mode.py /
build_protonation_state.py (no hardcoded atom index -- indices differ
between conformers/protonation states, e.g. after H removal upon
deprotonation).

Usage:
    python ~/work_turbomole/scripts/compute_b1_population_features.py \
        --from-csv layer_b/initial_labeled_subset.csv \
        --population-jobs layer_b/population_jobs \
        --xyz-dir layer_a/conformer_jobs --xyz-name structure.xyz \
        --output layer_b/population_feature.csv
"""

import argparse
import csv
import os
import re

import numpy as np

CARBONYL_BOND_MAX = 1.35
CO_BOND_MAX = 1.75
CH_BOND_MAX = 1.15

ATOMIC_NUMBER = {"h": 1, "c": 6, "n": 7, "o": 8}

SEN_RE = re.compile(
    r"shared electron number for the pair\s+(\d+)\s*([a-zA-Z]+)\s*-\s*(\d+)\s*([a-zA-Z]+)\s*=\s*([-\d.]+)"
)


def read_xyz(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0].split()[0])
    atoms = []
    for line in lines[2:2 + n]:
        parts = line.split()
        atoms.append((parts[0], np.array(list(map(float, parts[1:4])))))
    return atoms


def find_carbonyl_pair(atoms):
    """Same detection as identify_co_stretch_mode.py -- returns 0-based
    (o_idx, c_idx)."""
    def dist(i, j):
        return np.linalg.norm(atoms[i][1] - atoms[j][1])

    for oi, (sym_o, _) in enumerate(atoms):
        if sym_o != "O":
            continue
        if any(atoms[h][0] == "H" and dist(oi, h) < CH_BOND_MAX
               for h in range(len(atoms))):
            continue
        for ci, (sym_c, _) in enumerate(atoms):
            if sym_c != "C" or dist(oi, ci) >= CARBONYL_BOND_MAX:
                continue
            other_o = [j for j, (s, _) in enumerate(atoms)
                       if s == "O" and j != oi and dist(ci, j) < CO_BOND_MAX]
            if other_o:
                return oi, ci
    return None, None


def read_bond_order(pop_path, o_1based, c_1based):
    with open(pop_path) as f:
        for line in f:
            m = SEN_RE.search(line)
            if not m:
                continue
            idx1, _, idx2, _, value = m.groups()
            idx1, idx2 = int(idx1), int(idx2)
            if {idx1, idx2} == {o_1based, c_1based}:
                return float(value)
    return None


def read_bader_charges(acf_path, coord_path, o_1based, c_1based):
    symbols = []
    reading = False
    with open(coord_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("$coord"):
                reading = True
                continue
            if reading and (stripped.startswith("$end") or
                             (stripped.startswith("$") and stripped != "$coord")):
                break
            if reading and stripped:
                symbols.append(stripped.split()[3].lower())

    populations = []
    with open(acf_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                populations.append(float(parts[4]))

    def net_charge(idx_1based):
        sym = symbols[idx_1based - 1]
        z = ATOMIC_NUMBER.get(sym)
        pop = populations[idx_1based - 1]
        return z - pop if z is not None else float("nan")

    return net_charge(o_1based), net_charge(c_1based)


def main():
    parser = argparse.ArgumentParser(
        description="Extract C=O bond order (SEN) and Bader charges of the "
                    "carbonyl C/O atoms for a set of conformers.")
    parser.add_argument("--from-csv", required=True,
                         help="CSV with a 'job_dir' column")
    parser.add_argument("--population-jobs", required=True,
                         help="Parent folder of job_NNNN/ with pop.out + ACF.dat")
    parser.add_argument("--xyz-dir", required=True,
                         help="Parent folder of job_NNNN/ with the geometry "
                              "used to identify the carbonyl atoms")
    parser.add_argument("--xyz-name", default="structure.xyz")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.from_csv) as f:
        job_dirs = [row["job_dir"] for row in csv.DictReader(f)]

    results = []
    for jd in job_dirs:
        xyz_path = os.path.join(args.xyz_dir, jd, args.xyz_name)
        pop_dir = os.path.join(args.population_jobs, jd)
        pop_path = os.path.join(pop_dir, "pop.out")
        acf_path = os.path.join(pop_dir, "ACF.dat")
        coord_path = os.path.join(pop_dir, "coord")

        if not (os.path.isfile(pop_path) and os.path.isfile(acf_path)):
            print(f"WARNING: missing pop.out/ACF.dat for {jd}, skipping")
            continue

        atoms = read_xyz(xyz_path)
        o_idx, c_idx = find_carbonyl_pair(atoms)
        if o_idx is None:
            print(f"WARNING: no carbonyl found for {jd}, skipping")
            continue
        o_1b, c_1b = o_idx + 1, c_idx + 1

        bond_order = read_bond_order(pop_path, o_1b, c_1b)
        if bond_order is None:
            print(f"WARNING: SEN for O{o_1b}-C{c_1b} not found in {pop_path}, skipping")
            continue

        bader_o, bader_c = read_bader_charges(acf_path, coord_path, o_1b, c_1b)
        results.append((jd, bond_order, bader_o, bader_c))
        print(f"{jd}: SEN(O{o_1b}-C{c_1b})={bond_order:.4f}, "
              f"Bader(O{o_1b})={bader_o:+.4f}, Bader(C{c_1b})={bader_c:+.4f}")

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["job_dir", "co_bond_order_sen",
                          "bader_charge_carbonyl_o", "bader_charge_carbonyl_c"])
        for jd, bo, bo_charge, bc_charge in results:
            writer.writerow([jd, f"{bo:.6f}", f"{bo_charge:.6f}", f"{bc_charge:.6f}"])

    print(f"\nExtracted population features for {len(results)} / {len(job_dirs)} conformers.")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
