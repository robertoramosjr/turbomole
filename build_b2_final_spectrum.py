#!/usr/bin/env python3
"""
build_b2_final_spectrum.py

Builds the final B2 population-weighted UV-Vis envelope per
protonation state, from build_b2_consolidated_dataset.py's single
(E, f) per conformer (not a full multi-state escf.out table like
build_population_weighted_spectrum.py uses for Layer A's screening
spectra -- B2's per-conformer value IS the S1 prediction/label
itself). Same Gaussian-broadening convention as
build_population_weighted_spectrum.py (--broadening default 0.4 eV)
for consistency with the rest of the project.

Scope (protonation_and_conformer_investigation.md Sec. 9.28): only
neutral and monoanion_b are built here -- these are the two states
with either a validated per-conformer B2 model (monoanion_b) or a
corrected, calibrated surrogate + documented residual risk (neutral).
monoanion_a (only 26.6%% population coverage, no validated B2 model
for the rest) and dianion (N=3, out of scope for any model) are not
included in this final spectrum -- see build_b2_consolidated_dataset.py.

Usage:
    python ~/work_turbomole/scripts/build_b2_final_spectrum.py \
        --consolidated-dataset b2_consolidated_dataset.csv \
        --states neutral monoanion_b \
        --broadening 0.4 --npoints 3000 \
        --output-prefix b2_final_spectrum
"""

import argparse
import csv

import numpy as np


def gaussian_broaden(sticks, sigma, npoints, xmin, xmax):
    grid = np.linspace(xmin, xmax, npoints)
    profile = np.zeros_like(grid)
    for x0, height in sticks:
        profile += height * np.exp(-((grid - x0) ** 2) / (2 * sigma ** 2))
    return grid, profile


def main():
    parser = argparse.ArgumentParser(
        description="Build the population-weighted Gaussian-broadened UV-Vis "
                    "envelope per protonation state from the B2 consolidated "
                    "dataset.")
    parser.add_argument("--consolidated-dataset", required=True)
    parser.add_argument("--states", nargs="+", required=True)
    parser.add_argument("--broadening", type=float, default=0.4,
                         help="Gaussian sigma in eV [default: 0.4, matches "
                              "build_population_weighted_spectrum.py]")
    parser.add_argument("--npoints", type=int, default=3000)
    parser.add_argument("--xmin", type=float, default=None)
    parser.add_argument("--xmax", type=float, default=None)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    with open(args.consolidated_dataset) as f:
        all_rows = list(csv.DictReader(f))

    by_state = {}
    for r in all_rows:
        if r["protonation_state"] not in args.states:
            continue
        if not r["energy_eV"] or not r["osc_f"] or not r["population"]:
            continue
        by_state.setdefault(r["protonation_state"], []).append(
            (float(r["energy_eV"]), float(r["osc_f"]), float(r["population"])))

    all_x = [e for sticks in by_state.values() for e, _, _ in sticks]
    pad = args.broadening * 5
    xmin = args.xmin if args.xmin is not None else max(0.0, min(all_x) - pad)
    xmax = args.xmax if args.xmax is not None else max(all_x) + pad

    for state, sticks in by_state.items():
        total_pop = sum(p for _, _, p in sticks)
        weighted_sticks = [(e, f_osc * p) for e, f_osc, p in sticks]
        grid, profile = gaussian_broaden(weighted_sticks, args.broadening,
                                           args.npoints, xmin, xmax)
        out_path = f"{args.output_prefix}_{state}.dat"
        with open(out_path, "w") as f:
            f.write("descriptor energy_eV,absorption_weighted\n")
            for x, y in zip(grid, profile):
                f.write(f"{x:.6f} {y:.8f}\n")
        print(f"{state}: N={len(sticks)} conformers, population sum "
              f"(coverage)={total_pop:.4f}, peak absorption="
              f"{profile.max():.5f} at {grid[np.argmax(profile)]:.3f} eV "
              f"-> {out_path}")

    print(f"\nCommon grid: [{xmin:.3f}, {xmax:.3f}] eV, {args.npoints} points, "
          f"sigma={args.broadening} eV")


if __name__ == "__main__":
    main()
