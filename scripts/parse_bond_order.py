#!/usr/bin/env python3
"""
parse_bond_order.py

Parses "shared electron number" (SEN) pairs from a Turbomole
Mulliken population analysis output ($pop mulliken, printed via
dscf/ridft -proper), producing a Veusz-importable dataset.

SEN is a Mulliken-derived two-center bonding descriptor -- NOT
numerically equivalent to VASP+LOBSTER's COHP/ICOHP used in Januario &
Cabral (2026), but conceptually analogous (a measure of shared
electron density / bond strength between specific atom pairs). Used
here as the closest native Turbomole equivalent, since LOBSTER-style
projection is unnecessary for an already-localized (Gaussian) basis.

Usage:
    python ~/scripts/parse_bond_order.py --pop pop.out --output bond_order_dataset.dat
"""

import argparse
import re

ROW_RE = re.compile(
    r"shared electron number for the pair\s+(\d+)\s*([a-zA-Z]+)\s*-\s*(\d+)\s*([a-zA-Z]+)\s*=\s*([-\d.]+)"
)


def parse_sen(pop_path):
    pairs = []
    with open(pop_path, "r") as f:
        for line in f:
            match = ROW_RE.search(line)
            if match:
                idx1, elem1, idx2, elem2, value = match.groups()
                pairs.append({
                    "atom1_index": int(idx1),
                    "atom1_element": elem1,
                    "atom2_index": int(idx2),
                    "atom2_element": elem2,
                    "shared_electron_number": float(value),
                })

    if not pairs:
        raise ValueError(f"No 'shared electron number' pairs matched in {pop_path}.")

    return pairs


def main():
    parser = argparse.ArgumentParser(
        description="Parse Turbomole shared electron number (SEN) pairs into a Veusz dataset.")
    parser.add_argument("--pop", required=True, help="Path to pop.out (ridft/dscf -proper output)")
    parser.add_argument("--output", default="bond_order_dataset.dat",
                         help="Output .dat filename (Veusz descriptor format)")
    args = parser.parse_args()

    pairs = parse_sen(args.pop)
    print(f"Parsed {len(pairs)} atom-pair SEN values from {args.pop}")

    pairs_sorted = sorted(pairs, key=lambda p: p["shared_electron_number"], reverse=True)
    print("Top 5 strongest pairs:")
    for p in pairs_sorted[:5]:
        print(f"  {p['atom1_index']}{p['atom1_element']} - "
              f"{p['atom2_index']}{p['atom2_element']} = {p['shared_electron_number']:.4f}")

    with open(args.output, "w") as f:
        f.write("descriptor pair_label,shared_electron_number\n")
        for i, p in enumerate(pairs, start=1):
            label = f"{p['atom1_index']}{p['atom1_element']}-{p['atom2_index']}{p['atom2_element']}"
            f.write(f"{i} {p['shared_electron_number']:.6f}\n")

    # Also save a lookup CSV with the text labels, since Veusz text datasets
    # need separate handling from numeric ones for axis tick labels.
    with open(args.output.replace(".dat", "_labels.csv"), "w") as f:
        f.write("pair_index,pair_label,atom1_index,atom1_element,atom2_index,atom2_element,shared_electron_number\n")
        for i, p in enumerate(pairs, start=1):
            label = f"{p['atom1_index']}{p['atom1_element']}-{p['atom2_index']}{p['atom2_element']}"
            f.write(f"{i},{label},{p['atom1_index']},{p['atom1_element']},"
                     f"{p['atom2_index']},{p['atom2_element']},{p['shared_electron_number']:.6f}\n")

    print(f"Veusz dataset saved -> {args.output}")
    print(f"Label lookup saved -> {args.output.replace('.dat', '_labels.csv')}")


if __name__ == "__main__":
    main()