#!/usr/bin/env python3
"""
parse_becke_pdos.py

Converts a Multiwfn DOS_curve.txt export (Becke-partitioned total DOS and
element-resolved PDOS, already Gaussian-broadened) into a Veusz-importable
dataset. Unlike the original Mulliken-based PDOS (Sec. sec:dos), the Becke
partitioning scheme is real-space-based and yields strictly non-negative
element-resolved contributions.

Usage:
    python ~/scripts/parse_becke_pdos.py --curve DOS_curve.txt \
        --output dos_becke_dataset.dat
"""

import argparse

HARTREE_TO_EV = 27.211386245988


def parse_dos_curve(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            energy_au, tdos, pdos_o, pdos_c, pdos_h = (
                float(parts[0]), float(parts[1]), float(parts[2]),
                float(parts[3]), float(parts[4])
            )
            rows.append((energy_au, tdos, pdos_o, pdos_c, pdos_h))
    if not rows:
        raise ValueError(f"No data rows parsed from {path}.")
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Convert a Multiwfn Becke-PDOS export to a Veusz dataset.")
    parser.add_argument("--curve", required=True, help="Path to DOS_curve.txt")
    parser.add_argument("--output", default="dos_becke_dataset.dat")
    args = parser.parse_args()

    rows = parse_dos_curve(args.curve)
    print(f"Parsed {len(rows)} points from {args.curve}")

    negatives = [r for r in rows if min(r[1:]) < 0]
    if negatives:
        print(f"WARNING: {len(negatives)} points still have a negative "
              f"component -- check the Multiwfn export before trusting this dataset.")
    else:
        print("Confirmed: all TDOS and PDOS values are non-negative.")

    with open(args.output, "w") as f:
        f.write("descriptor energy_eV,TDOS,pDOS_O,pDOS_C,pDOS_H\n")
        for energy_au, tdos, pdos_o, pdos_c, pdos_h in rows:
            energy_ev = energy_au * HARTREE_TO_EV
            f.write(f"{energy_ev:.6f} {tdos:.6f} {pdos_o:.6f} "
                     f"{pdos_c:.6f} {pdos_h:.6f}\n")

    print(f"Veusz dataset saved -> {args.output}")


if __name__ == "__main__":
    main()