#!/usr/bin/env python3
"""
compute_b1_dihedral_feature.py

Computes the carbonyl-vinyl conjugation dihedral
O(carbonyl)=C(carbonyl)-C(alpha)=C(beta) for every conformer in a
job-directory table -- a second candidate feature for B1's GPR,
alongside the O...H distance (compute_b1_features.py).

Rationale: artepilin C's carboxyl is the end of an acrylic-acid-type
arm (Ar-CH=CH-COOH), confirmed by geometry (C-alpha to carbonyl-C:
~1.465 A single bond; C-alpha to C-beta: ~1.353 A double bond, on the
Paper 1 reference conformer). How coplanar this arm is with the C=O
controls pi-conjugation between the vinyl double bond and the
carbonyl -- more planar means stronger conjugation, which delocalizes
electron density out of the C=O bond and should shift its stretching
frequency; a twisted (non-planar) arrangement breaks that conjugation.
This dihedral is a natural second feature to test given the O...H
distance alone didn't carry consistent signal across all 4 protonation
states (see train_b1_gpr.py's leave-one-out results).

Atom identification is entirely geometry-based (same bond-distance
philosophy as the rest of this project's scripts, not hardcoded
indices): finds the carbonyl O/C pair the same way as
identify_co_stretch_mode.py, then the carbonyl carbon's third
substituent that is itself a carbon (C-alpha), then C-alpha's shortest
C-C bond excluding the carbonyl carbon (C-beta, the double-bond
partner). This assumes the specific connectivity confirmed for
artepilin C (an alpha,beta-unsaturated carboxyl) -- verify with
--list-atoms on a new molecule before trusting this blindly.

Usage:
    python ~/work_turbomole/scripts/compute_b1_dihedral_feature.py \
        --from-csv layer_b/initial_labeled_subset.csv \
        --xyz-dir layer_a/conformer_jobs --xyz-name structure.xyz \
        --output layer_b/dihedral_feature.csv
"""

import argparse
import csv
import os

import numpy as np

CARBONYL_BOND_MAX = 1.35
CO_BOND_MAX = 1.75
CH_BOND_MAX = 1.15
CC_BOND_MAX = 1.65


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


def find_conjugation_chain(atoms, o_idx, c_idx):
    """Return (c_alpha, c_beta) 0-based indices for the vinyl arm
    attached to the carbonyl carbon."""
    def dist(i, j):
        return np.linalg.norm(atoms[i][1] - atoms[j][1])

    c_alpha = None
    for j, (s, _) in enumerate(atoms):
        if s == "C" and j != c_idx and dist(c_idx, j) < CC_BOND_MAX:
            c_alpha = j
            break
    if c_alpha is None:
        return None, None

    c_beta = None
    best_d = None
    for j, (s, _) in enumerate(atoms):
        if s != "C" or j in (c_idx, c_alpha):
            continue
        d = dist(c_alpha, j)
        if d < CC_BOND_MAX and (best_d is None or d < best_d):
            c_beta, best_d = j, d
    return c_alpha, c_beta


def dihedral_angle(p0, p1, p2, p3):
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def main():
    parser = argparse.ArgumentParser(
        description="Compute the O=C-C=C conjugation dihedral for the "
                    "carbonyl-vinyl arm.")
    parser.add_argument("--from-csv", required=True)
    parser.add_argument("--xyz-dir", required=True)
    parser.add_argument("--xyz-name", default="structure.xyz")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.from_csv) as f:
        job_dirs = [row["job_dir"] for row in csv.DictReader(f)]

    results = []
    for jd in job_dirs:
        xyz_path = os.path.join(args.xyz_dir, jd, args.xyz_name)
        atoms = read_xyz(xyz_path)
        o_idx, c_idx = find_carbonyl_pair(atoms)
        if o_idx is None:
            print(f"WARNING: no carbonyl found for {jd}")
            continue
        c_alpha, c_beta = find_conjugation_chain(atoms, o_idx, c_idx)
        if c_alpha is None or c_beta is None:
            print(f"WARNING: no conjugated vinyl arm found for {jd}")
            continue
        angle = dihedral_angle(atoms[o_idx][1], atoms[c_idx][1],
                                atoms[c_alpha][1], atoms[c_beta][1])
        # fold to [0, 90] -- planarity (0 or 180) is what matters for
        # conjugation strength, not the arbitrary sign/direction of twist
        planarity_deviation = min(abs(angle), abs(180 - abs(angle)))
        results.append((jd, angle, planarity_deviation))

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["job_dir", "oc_dihedral_deg", "planarity_deviation_deg"])
        for jd, a, p in results:
            writer.writerow([jd, f"{a:.2f}", f"{p:.2f}"])

    devs = [p for _, _, p in results]
    print(f"Computed conjugation dihedral for {len(results)} / {len(job_dirs)} conformers.")
    print(f"Planarity deviation range: {min(devs):.2f} - {max(devs):.2f} deg "
          f"(0 = fully planar/conjugated, 90 = fully twisted/decoupled)")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
