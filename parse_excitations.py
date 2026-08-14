#!/usr/bin/env python3
"""
parse_excitations.py

Parses the "SUMMARY OF EXCITATION ENERGIES AND DIPOLE OSCILLATOR
STRENGTHS" table from a Turbomole escf.out file (TD-DFT/RPA singlet
excitations), producing:

  1. A Veusz-importable "stick spectrum" dataset (excitation energy vs.
     oscillator strength) -- raw data behind the optical absorption
     spectrum (Fig. optical).
  2. An oscillator-strength-weighted broadened envelope
     (jdos_dataset.dat) -- used ONLY in the optical absorption figure,
     not the DOS figure.
  3. An equal-weight excitation-energy density (each TD-DFT root
     contributes with unit weight, regardless of oscillator strength)
     -- used in the DOS figure panel (b), per the paper's current
     wording ("equal-weight TD-DFT excitation-energy density").

Usage:
    python ~/scripts/parse_excitations.py --escf escf.out \
        --output-sticks excitation_sticks.dat \
        --output-jdos jdos_dataset.dat \
        --output-density excitation_density_dataset.dat \
        --broadening 0.4 --npoints 2000
"""

import argparse
import re
import numpy as np

TABLE_ROW_RE = re.compile(
    r"^\s*\|\s*(\d+)\s+a\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
)


def parse_escf_table(escf_path):
    rows = []
    with open(escf_path, "r") as f:
        for line in f:
            match = TABLE_ROW_RE.match(line)
            if match:
                state, e_eh, e_ev, e_cm1, e_nm, osc_vel, osc_len = match.groups()
                rows.append({
                    "state": int(state),
                    "energy_Eh": float(e_eh),
                    "energy_eV": float(e_ev),
                    "energy_cm1": float(e_cm1),
                    "energy_nm": float(e_nm),
                    "osc_vel": float(osc_vel),
                    "osc_len": float(osc_len),
                })

    if not rows:
        raise ValueError(
            f"No excitation table rows matched in {escf_path}. "
            f"Check that the file contains the 'SUMMARY OF EXCITATION ENERGIES' "
            f"table and that the row format matches what this script expects."
        )

    return rows


def write_stick_spectrum(rows, output_path):
    with open(output_path, "w") as f:
        f.write("descriptor state,energy_eV,osc_len,osc_vel\n")
        for r in rows:
            f.write(f"{r['state']} {r['energy_eV']:.6f} "
                     f"{r['osc_len']:.6f} {r['osc_vel']:.6f}\n")
    print(f"Stick spectrum saved -> {output_path}")


def gaussian_broaden_weighted(rows, broadening_eV, npoints, emin=None, emax=None):
    """Oscillator-strength-weighted envelope -- for the optical absorption figure only."""
    energies = np.array([r["energy_eV"] for r in rows])
    oscillators = np.array([r["osc_len"] for r in rows])

    if emin is None:
        emin = max(0.0, energies.min() - 5 * broadening_eV)
    if emax is None:
        emax = energies.max() + 5 * broadening_eV

    grid = np.linspace(emin, emax, npoints)
    profile = np.zeros_like(grid)
    for e0, f in zip(energies, oscillators):
        profile += f * np.exp(-((grid - e0) ** 2) / (2 * broadening_eV ** 2))

    return grid, profile


def gaussian_broaden_equal_weight(rows, broadening_eV, npoints, emin=None, emax=None):
    """Equal-weight excitation-energy density -- for the DOS figure panel (b)."""
    energies = np.array([r["energy_eV"] for r in rows])

    if emin is None:
        emin = max(0.0, energies.min() - 5 * broadening_eV)
    if emax is None:
        emax = energies.max() + 5 * broadening_eV

    grid = np.linspace(emin, emax, npoints)
    profile = np.zeros_like(grid)
    for e0 in energies:
        profile += np.exp(-((grid - e0) ** 2) / (2 * broadening_eV ** 2))

    return grid, profile


def write_jdos_dataset(grid, profile, output_path):
    with open(output_path, "w") as f:
        f.write("descriptor jdos_energy_eV,jdos_intensity\n")
        for e, i in zip(grid, profile):
            f.write(f"{e:.6f} {i:.6f}\n")
    print(f"Oscillator-strength-weighted envelope saved -> {output_path}")


def write_density_dataset(grid, profile, output_path):
    with open(output_path, "w") as f:
        f.write("descriptor density_energy_eV,excitation_density\n")
        for e, i in zip(grid, profile):
            f.write(f"{e:.6f} {i:.6f}\n")
    print(f"Equal-weight excitation-energy density saved -> {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Parse Turbomole escf.out excitation table into Veusz-importable datasets.")
    parser.add_argument("--escf", required=True, help="Path to escf.out")
    parser.add_argument("--output-sticks", default="excitation_sticks.dat")
    parser.add_argument("--output-jdos", default="jdos_dataset.dat",
                         help="Oscillator-strength-weighted envelope (optical figure)")
    parser.add_argument("--output-density", default="excitation_density_dataset.dat",
                         help="Equal-weight excitation-energy density (DOS figure)")
    parser.add_argument("--broadening", type=float, default=0.4)
    parser.add_argument("--npoints", type=int, default=2000)
    args = parser.parse_args()

    rows = parse_escf_table(args.escf)
    print(f"Parsed {len(rows)} excited states "
          f"(lowest: {rows[0]['energy_eV']:.3f} eV, "
          f"highest: {rows[-1]['energy_eV']:.3f} eV)")

    write_stick_spectrum(rows, args.output_sticks)

    grid_w, profile_w = gaussian_broaden_weighted(rows, args.broadening, args.npoints)
    write_jdos_dataset(grid_w, profile_w, args.output_jdos)

    grid_eq, profile_eq = gaussian_broaden_equal_weight(rows, args.broadening, args.npoints)
    write_density_dataset(grid_eq, profile_eq, args.output_density)


if __name__ == "__main__":
    main()