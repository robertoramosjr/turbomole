#!/usr/bin/env python3
"""
parse_qpenergies_dos.py

Parses a Turbomole qpenergies.dat file (G0W0 quasiparticle corrections)
and builds a Gaussian-broadened "orbital density of states" for both the
Kohn-Sham (HSE06) and GW-corrected eigenvalue sets, on the same energy
grid, for direct overlay -- e.g. Fig. 3 (HSE06 DOS) vs. the QP-corrected
spectrum requested by the coauthor.

This is deliberately an UNPROJECTED (element-agnostic) orbital density --
each orbital contributes with unit weight, the same equal-weight
broadening convention already used for the TD-DFT excitation density
(parse_excitations.py, --output-density). It sidesteps the Mulliken
population-weighting issue entirely, since no atom projection is
involved here.

Usage:
    python ~/scripts/parse_qpenergies_dos.py --qpenergies qpenergies.dat \
        --output-ks dos_ks_dataset.dat --output-qp dos_qp_dataset.dat \
        --broadening 0.136 --npoints 4000 --emin -20 --emax 10
"""

import argparse
import re
import numpy as np

ROW_RE = re.compile(
    r"^\s*(\d+)([a-zA-Z]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$"
)


def parse_qpenergies(path):
    """Parse orbital index, KS eigenvalue, and G0W0 QP eigenvalue (eV)."""
    orbitals = []
    with open(path, "r") as f:
        for line in f:
            if line.strip().startswith("$"):
                continue
            match = ROW_RE.match(line)
            if match:
                index, irrep, e_ks, e_qp, _col3, _sigma_x, _sigma_c = match.groups()
                orbitals.append({
                    "index": int(index),
                    "irrep": irrep,
                    "e_ks_eV": float(e_ks),
                    "e_qp_eV": float(e_qp),
                })

    if not orbitals:
        raise ValueError(
            f"No orbital rows matched in {path}. Check that the file contains "
            f"the '$qpenergies ... Kohn-Sham E. ... G0W0 ... QPs' table and "
            f"that the row format matches what this script expects."
        )

    return orbitals


def gaussian_broadened_density(energies_eV, broadening_eV, npoints, emin, emax):
    """Equal-weight orbital density: each orbital contributes one unit-area Gaussian."""
    grid = np.linspace(emin, emax, npoints)
    profile = np.zeros_like(grid)
    for e0 in energies_eV:
        profile += np.exp(-((grid - e0) ** 2) / (2 * broadening_eV ** 2))
    return grid, profile


def write_density_dataset(grid, profile, label, output_path):
    with open(output_path, "w") as f:
        f.write(f"descriptor energy_eV,{label}_density\n")
        for e, i in zip(grid, profile):
            f.write(f"{e:.6f} {i:.6f}\n")
    print(f"{label} orbital density saved -> {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Build Gaussian-broadened KS vs. G0W0-corrected orbital "
                    "density datasets from a Turbomole qpenergies.dat file, "
                    "for direct overlay in Veusz.")
    parser.add_argument("--qpenergies", required=True, help="Path to qpenergies.dat")
    parser.add_argument("--output-ks", default="dos_ks_dataset.dat",
                         help="Output dataset for the Kohn-Sham (HSE06) orbital density")
    parser.add_argument("--output-qp", default="dos_qp_dataset.dat",
                         help="Output dataset for the G0W0-corrected orbital density")
    parser.add_argument("--broadening", type=float, default=0.136,
                         help="Gaussian FWHM in eV (default 0.136 eV = 0.005 a.u., "
                              "matching the tightened DOS width already validated "
                              "for Fig. 3)")
    parser.add_argument("--npoints", type=int, default=4000)
    parser.add_argument("--emin", type=float, default=-20.0)
    parser.add_argument("--emax", type=float, default=10.0)
    args = parser.parse_args()

    orbitals = parse_qpenergies(args.qpenergies)
    print(f"Parsed {len(orbitals)} orbitals from {args.qpenergies}")

    n_in_window = sum(1 for o in orbitals if args.emin <= o["e_ks_eV"] <= args.emax)
    print(f"{n_in_window} orbitals fall inside the plotted window "
          f"[{args.emin}, {args.emax}] eV (Kohn-Sham energy)")
    if n_in_window == 0:
        print("WARNING: no orbitals fall inside the plotted window. Check units "
              "(this script assumes eV) and the --emin/--emax range.")

    e_ks = [o["e_ks_eV"] for o in orbitals]
    e_qp = [o["e_qp_eV"] for o in orbitals]

    grid_ks, profile_ks = gaussian_broadened_density(
        e_ks, args.broadening, args.npoints, args.emin, args.emax)
    write_density_dataset(grid_ks, profile_ks, "ks", args.output_ks)

    grid_qp, profile_qp = gaussian_broadened_density(
        e_qp, args.broadening, args.npoints, args.emin, args.emax)
    write_density_dataset(grid_qp, profile_qp, "qp", args.output_qp)


if __name__ == "__main__":
    main()