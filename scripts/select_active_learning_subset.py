#!/usr/bin/env python3
"""
select_active_learning_subset.py

Chooses the initial labeled subset for Layer B's active-learning
bootstrap (B1/B2 GPR surrogates), out of Layer A's population-filtered
conformers of one protonation state.

Selection combines two criteria, not one:
  1. Representativeness -- the top --n-population-seed conformers by
     population are always included. These matter most for the final
     Boltzmann-weighted properties, so the surrogate should never be
     extrapolating for them.
  2. Diversity -- the remaining slots (up to --n-target total) are
     filled by farthest-point sampling (FPS) on structural RMSD
     (Kabsch-aligned Cartesian coordinates), seeded from the
     population picks: repeatedly add whichever remaining conformer
     has the largest RMSD to its nearest already-selected neighbor.
     This is what active learning needs to bootstrap well -- a GPR
     trained only on the lowest-energy cluster has no information
     about how the delta-learning correction behaves elsewhere in
     conformational space, and its uncertainty estimates (used to pick
     the next active-learning labels) would be unreliable outside that
     cluster.

If the whole population-filtered ensemble is already <= --n-target,
everything is selected (no algorithm needed) -- this is expected for
the smaller protonation states (e.g. dianion had only 3 survivors in
this project).

Usage:
    python ~/work_turbomole/scripts/select_active_learning_subset.py \
        --population-table layer_a/population_table_dft.csv \
        --conformer-jobs layer_a/conformer_jobs \
        --n-target 20 --n-population-seed 3 \
        --output layer_b/initial_labeled_subset.csv
"""

import argparse
import csv
import os

import numpy as np


def read_coord(path):
    """Read a Turbomole 'coord' file, return (symbols, xyz in Bohr)."""
    symbols = []
    coords = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 4 and parts[0].replace(".", "").replace(
                    "-", "").isdigit():
                x, y, z, sym = parts
                coords.append([float(x), float(y), float(z)])
                symbols.append(sym.lower())
    return symbols, np.array(coords)


def kabsch_rmsd(p, q):
    """RMSD between two (N,3) coordinate sets after optimal
    superposition (Kabsch algorithm) -- assumes matching atom order,
    which holds here since every conformer of a given protonation
    state comes from the same starting topology."""
    p = p - p.mean(axis=0)
    q = q - q.mean(axis=0)
    h = p.T @ q
    u, s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    corr = np.diag([1, 1, d])
    r = vt.T @ corr @ u.T
    p_aligned = (r @ p.T).T
    return np.sqrt(np.mean(np.sum((p_aligned - q) ** 2, axis=1)))


def main():
    parser = argparse.ArgumentParser(
        description="Select the initial labeled subset for Layer B "
                    "active learning: top-population conformers plus "
                    "farthest-point structural diversity sampling.")
    parser.add_argument("--population-table", required=True,
                         help="Layer A's population_table_dft.csv")
    parser.add_argument("--conformer-jobs", required=True,
                         help="Parent folder of job_NNNN/coord dirs "
                              "(r2SCAN-3c-optimized geometries)")
    parser.add_argument("--n-target", type=int, default=20,
                         help="Total labeled subset size [default: 20, "
                              "within the 15-30 range decided for this project]")
    parser.add_argument("--n-population-seed", type=int, default=3,
                         help="How many top-population conformers are always "
                              "included before diversity sampling [default: 3]")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.population_table) as f:
        # population_table_dft.csv (rerank_layer_a_dft.py) already
        # contains only Layer A's survivors -- no "kept" column to filter on.
        rows = list(csv.DictReader(f))
    rows.sort(key=lambda r: -float(r["r2scan3c_population"]))
    job_dirs = [r["job_dir"] for r in rows]
    populations = {r["job_dir"]: float(r["r2scan3c_population"]) for r in rows}

    n_available = len(job_dirs)
    if n_available <= args.n_target:
        selected = list(job_dirs)
        reasons = {jd: "all_survivors_below_target" for jd in selected}
        print(f"Only {n_available} survivors (<= target {args.n_target}) -- "
              f"selecting all of them, no sampling needed.")
    else:
        coords = {
            jd: read_coord(os.path.join(args.conformer_jobs, jd, "coord"))[1]
            for jd in job_dirs
        }

        n_seed = min(args.n_population_seed, args.n_target, n_available)
        selected = job_dirs[:n_seed]
        reasons = {jd: "population_seed" for jd in selected}
        remaining = [jd for jd in job_dirs if jd not in selected]

        while len(selected) < args.n_target and remaining:
            best_jd, best_min_rmsd = None, -1.0
            for jd in remaining:
                min_rmsd = min(
                    kabsch_rmsd(coords[jd], coords[s]) for s in selected
                )
                if min_rmsd > best_min_rmsd:
                    best_jd, best_min_rmsd = jd, min_rmsd
            selected.append(best_jd)
            reasons[best_jd] = f"diversity_fps(rmsd={best_min_rmsd:.3f}A)"
            remaining.remove(best_jd)

        print(f"Selected {len(selected)} / {n_available} conformers: "
              f"{n_seed} by population, {len(selected) - n_seed} by "
              f"structural diversity (farthest-point sampling).")

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["job_dir", "population", "selection_reason"])
        for jd in selected:
            writer.writerow([jd, f"{populations[jd]:.6f}", reasons[jd]])

    total_pop_covered = sum(populations[jd] for jd in selected)
    print(f"Total population covered by this labeled subset: "
          f"{total_pop_covered:.1%}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
