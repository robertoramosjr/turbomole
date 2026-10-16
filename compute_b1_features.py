#!/usr/bin/env python3
"""
compute_b1_features.py

Computes the structural feature Paper 1 flagged as correlating with
the carbonyl band shift -- the O(carbonyl)...H(nearest donor) distance
(2.281 A in the Paper 1 reference conformer) -- for every conformer
listed in a job-directory table, as a candidate delta-learning feature
for B1's GPR. Carbonyl oxygen is detected the same geometry-based way
as build_protonation_state.py / identify_co_stretch_mode.py (no
hardcoded atom index).

"Nearest H donor" is simply the closest H atom to the carbonyl oxygen,
excluding H atoms bonded to the carbonyl's own carbon (there normally
are none in a carboxyl/carboxylate, but this guards against
mismeasuring a covalent bond as a hydrogen-bond donor distance).

Usage:
    python ~/work_turbomole/scripts/compute_b1_features.py \
        --job-list job_dirs.txt --xyz-name structure.xyz \
        --output layer_b/oh_distance_feature.csv

    # or directly against a labeled-subset/population CSV's job_dir column:
    python ~/work_turbomole/scripts/compute_b1_features.py \
        --from-csv layer_b/initial_labeled_subset.csv \
        --xyz-dir layer_a/conformer_jobs --xyz-name structure.xyz \
        --output layer_b/oh_distance_feature.csv
"""

import argparse
import csv
import os

import numpy as np

CARBONYL_BOND_MAX = 1.35
CO_BOND_MAX = 1.75
CH_BOND_MAX = 1.15


def read_xyz(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0].split()[0])
    atoms = []
    for line in lines[2:2 + n]:
        parts = line.split()
        atoms.append((parts[0], np.array(list(map(float, parts[1:4])))))
    return atoms


def find_carbonyl_o(atoms):
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


def nearest_h_distance(atoms, o_idx, c_idx):
    best = None
    for hi, (sym, pos) in enumerate(atoms):
        if sym != "H":
            continue
        if np.linalg.norm(atoms[c_idx][1] - pos) < CH_BOND_MAX:
            continue  # skip H covalently bonded to the carbonyl carbon
        d = np.linalg.norm(atoms[o_idx][1] - pos)
        if best is None or d < best:
            best = d
    return best


def main():
    parser = argparse.ArgumentParser(
        description="Compute the O(carbonyl)...H(nearest donor) distance "
                    "feature for a set of conformers.")
    parser.add_argument("--from-csv", required=True,
                         help="CSV with a 'job_dir' column (e.g. "
                              "initial_labeled_subset.csv or population_table_dft.csv)")
    parser.add_argument("--xyz-dir", required=True,
                         help="Parent folder of job_NNNN/ dirs")
    parser.add_argument("--xyz-name", default="structure.xyz")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.from_csv) as f:
        job_dirs = [row["job_dir"] for row in csv.DictReader(f)]

    results = []
    for jd in job_dirs:
        xyz_path = os.path.join(args.xyz_dir, jd, args.xyz_name)
        atoms = read_xyz(xyz_path)
        o_idx, c_idx = find_carbonyl_o(atoms)
        if o_idx is None:
            print(f"WARNING: no carbonyl found for {jd}")
            continue
        d = nearest_h_distance(atoms, o_idx, c_idx)
        results.append((jd, d))

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["job_dir", "oh_distance_angstrom"])
        for jd, d in results:
            writer.writerow([jd, f"{d:.4f}"])

    dists = [d for _, d in results]
    print(f"Computed O...H distance for {len(results)} / {len(job_dirs)} conformers.")
    print(f"Range: {min(dists):.3f} - {max(dists):.3f} A "
          f"(Paper 1 reference conformer: 2.281 A)")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
