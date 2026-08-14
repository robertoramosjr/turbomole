#!/usr/bin/env python3
"""
convergence_table.py

Builds a numerical-convergence comparison table (total SCF energy +
dominant TD-DFT excitation) across multiple runs that differ in one
numerical setting (integration grid, SCF threshold, basis set, etc.),
with the conformer/geometry held fixed. Reusable for any of the
"robustness" tests in item 4 (grid, SCF, basis) -- not just grid size.

Usage:
    python ~/scripts/convergence_table.py \
        --run "m3 (default)":artepillin_C_d3_gridm3/ridft_grid.out:artepillin_C_d3_gridm3/escf_grid.out \
        --run "m4 (production)":artepillin_C_d3/ridft_d3.out:artepillin_C_d3/escf_d3.out \
        --run "m5 (fine)":artepillin_C_d3_gridm5/ridft_grid.out:artepillin_C_d3_gridm5/escf_grid.out \
        --output grid_convergence.csv
"""

import argparse
import re

ENERGY_RE = re.compile(r"\|\s*total energy\s*=\s*([-\d.]+)\s*\|")
TABLE_ROW_RE = re.compile(
    r"^\s*\|\s*(\d+)\s+a\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
)


def parse_total_energy(ridft_path):
    """Returns the LAST 'total energy' value found (final converged SCF)."""
    energy = None
    with open(ridft_path, "r") as f:
        for line in f:
            match = ENERGY_RE.search(line)
            if match:
                energy = float(match.group(1))
    if energy is None:
        raise ValueError(f"No 'total energy' line found in {ridft_path}")
    return energy


def parse_dominant_excitation(escf_path):
    rows = []
    with open(escf_path, "r") as f:
        for line in f:
            match = TABLE_ROW_RE.match(line)
            if match:
                state, e_eh, e_ev, e_cm1, e_nm, osc_vel, osc_len = match.groups()
                rows.append({"state": int(state), "energy_eV": float(e_ev), "osc_len": float(osc_len)})
    if not rows:
        raise ValueError(f"No excitation rows found in {escf_path}")
    return max(rows, key=lambda r: r["osc_len"])


def main():
    parser = argparse.ArgumentParser(description="Build a numerical-convergence table (SCF energy + dominant excitation).")
    parser.add_argument("--run", action="append", required=True,
                        metavar="LABEL|RIDFT_OUT|ESCF_OUT",
                        help="Label, path to ridft.out, path to escf.out (pipe-separated). Repeat.")
    parser.add_argument("--output", default="convergence_table.csv")
    args = parser.parse_args()

    results = []
    for entry in args.run:
        label, ridft_path, escf_path = entry.split("|", 2)
        energy = parse_total_energy(ridft_path)
        dominant = parse_dominant_excitation(escf_path)
        results.append({"label": label, "total_energy_Eh": energy, **dominant})

    reference_energy = results[0]["total_energy_Eh"]
    reference_excitation = results[0]["energy_eV"]

    print(f"{'Run':<20} {'E_total (Eh)':<16} {'dE vs ref (Eh)':<16} "
          f"{'Dominant (eV)':<16} {'d(exc) vs ref (eV)':<18}")
    for r in results:
        d_energy = r["total_energy_Eh"] - reference_energy
        d_exc = r["energy_eV"] - reference_excitation
        print(f"{r['label']:<20} {r['total_energy_Eh']:<16.8f} {d_energy:<+16.8f} "
              f"{r['energy_eV']:<16.4f} {d_exc:<+18.4f}")

    with open(args.output, "w") as f:
        f.write("label,total_energy_Eh,dE_vs_reference_Eh,dominant_excitation_eV,d_excitation_vs_reference_eV\n")
        for r in results:
            d_energy = r["total_energy_Eh"] - reference_energy
            d_exc = r["energy_eV"] - reference_excitation
            f.write(f"{r['label']},{r['total_energy_Eh']:.8f},{d_energy:+.8f},"
                     f"{r['energy_eV']:.4f},{d_exc:+.4f}\n")

    print(f"\nTable saved -> {args.output}")


if __name__ == "__main__":
    main()