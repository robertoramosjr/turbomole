#!/usr/bin/env python3
"""
compute_rdf.py (v2 -- Gaussian-broadened)

Computes a Gaussian-broadened radial distribution function for
user-defined functional groups, referenced to each group's center of
mass. Replaces the raw fine-binned histogram (which produced sparse,
needle-like spikes dominated by near-r=0 normalization artifacts) with
a continuous kernel density estimate -- a standard technique for
producing smooth, presentable RDF curves from a small number of
pairwise distances.

Usage:
    python ~/scripts/compute_rdf.py --coord coord --groups groups.json \
        --sigma 0.15 --rmax 14.0 --npoints 500 --output rdf_dataset.dat
"""

import argparse
import json
import numpy as np

BOHR_TO_ANGSTROM = 0.529177210903

ATOMIC_MASS = {
    "h": 1.008, "c": 12.011, "n": 14.007, "o": 15.999,
    "s": 32.06, "cl": 35.45, "f": 18.998, "p": 30.974,
}


def read_turbomole_coord(coord_path):
    symbols = []
    coords_bohr = []
    reading = False
    with open(coord_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("$coord"):
                reading = True
                continue
            if reading and (stripped.startswith("$end") or
                             (stripped.startswith("$") and stripped != "$coord")):
                break
            if reading and stripped:
                x, y, z, symbol = stripped.split()[:4]
                coords_bohr.append([float(x), float(y), float(z)])
                symbols.append(symbol.lower())
    return symbols, np.array(coords_bohr) * BOHR_TO_ANGSTROM


def center_of_mass(symbols, coords, indices):
    masses = np.array([ATOMIC_MASS[symbols[i - 1]] for i in indices])
    positions = np.array([coords[i - 1] for i in indices])
    return np.average(positions, axis=0, weights=masses)


def gaussian_broadened_rdf(coords, reference_point, sigma, rmax, npoints):
    """
    Kernel density estimate: each pairwise distance contributes a
    Gaussian of width sigma, normalized by shell volume (4*pi*r^2) to
    preserve the standard g(r)-style radial density interpretation,
    without the small-bin-count artifact of raw histogramming.
    """
    distances = np.linalg.norm(coords - reference_point, axis=1)
    distances = distances[distances > 1e-6]  # exclude the reference atom itself, if present

    r_grid = np.linspace(0.01, rmax, npoints)
    density = np.zeros_like(r_grid)

    for d in distances:
        density += np.exp(-((r_grid - d) ** 2) / (2 * sigma ** 2))

    shell_volume = 4.0 * np.pi * r_grid ** 2
    g_r = density / shell_volume

    return r_grid, g_r


def main():
    parser = argparse.ArgumentParser(
        description="Compute Gaussian-broadened RDF histograms from a Turbomole coord file.")
    parser.add_argument("--coord", required=True)
    parser.add_argument("--groups", required=True)
    parser.add_argument("--sigma", type=float, default=0.15,
                         help="Gaussian broadening width (Angstrom)")
    parser.add_argument("--rmax", type=float, default=14.0)
    parser.add_argument("--npoints", type=int, default=500)
    parser.add_argument("--output", default="rdf_dataset.dat")
    args = parser.parse_args()

    symbols, coords = read_turbomole_coord(args.coord)

    with open(args.groups, "r") as f:
        groups = json.load(f)

    group_names = list(groups.keys())
    all_gr = {}
    reference_r = None

    for group_name, indices in groups.items():
        com = center_of_mass(symbols, coords, indices)
        r_grid, g_r = gaussian_broadened_rdf(coords, com, args.sigma, args.rmax, args.npoints)
        all_gr[group_name] = g_r
        reference_r = r_grid
        print(f"[{group_name}] center of mass (Angstrom): {com}")

    descriptor_names = ["r_angstrom"] + [f"gr_{name}" for name in group_names]
    with open(args.output, "w") as f:
        f.write("descriptor " + ",".join(descriptor_names) + "\n")
        for i, r in enumerate(reference_r):
            row = [f"{r:.4f}"] + [f"{all_gr[name][i]:.6f}" for name in group_names]
            f.write(" ".join(row) + "\n")

    print(f"Veusz dataset saved -> {args.output}")


if __name__ == "__main__":
    main()