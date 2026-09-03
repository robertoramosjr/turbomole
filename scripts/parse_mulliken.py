#!/usr/bin/env python3
"""
parse_mulliken.py

Parses the Mulliken atomic charges table from a Turbomole population
analysis output ($pop mulliken, via ridft/dscf -proper), producing a
Veusz-importable dataset. Intended as a cross-check against the Bader
net charges computed earlier (parse_bader.py / bader_dataset.dat).

Usage:
    python ~/scripts/parse_mulliken.py --pop pop.out --output mulliken_dataset.dat
"""

import argparse
import re

ROW_RE = re.compile(
    r"^\s*(\d+)([a-zA-Z]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
)


def parse_mulliken(pop_path):
    charges = []
    in_table = False

    with open(pop_path, "r") as f:
        for line in f:
            if "atomic populations from total density" in line:
                in_table = True
                continue
            if in_table:
                match = ROW_RE.match(line)
                if match:
                    idx, element, charge = match.group(1), match.group(2), match.group(3)
                    charges.append({
                        "atom_index": int(idx),
                        "element": element,
                        "mulliken_charge": float(charge),
                    })
                elif charges and line.strip() == "":
                    # blank line after the table -- stop capturing
                    break

    if not charges:
        raise ValueError(f"No Mulliken charge rows matched in {pop_path}.")

    return charges


def main():
    parser = argparse.ArgumentParser(
        description="Parse Turbomole Mulliken atomic charges into a Veusz dataset.")
    parser.add_argument("--pop", required=True, help="Path to pop.out")
    parser.add_argument("--output", default="mulliken_dataset.dat",
                         help="Output .dat filename (Veusz descriptor format)")
    args = parser.parse_args()

    charges = parse_mulliken(args.pop)
    print(f"Parsed {len(charges)} Mulliken atomic charges from {args.pop}")

    total = sum(c["mulliken_charge"] for c in charges)
    print(f"Sum of Mulliken charges (should be close to 0 for a neutral molecule): {total:.4f}")

    with open(args.output, "w") as f:
        f.write("descriptor atom_index,mulliken_charge\n")
        for c in charges:
            f.write(f"{c['atom_index']} {c['mulliken_charge']:.6f}\n")

    print(f"Veusz dataset saved -> {args.output}")


if __name__ == "__main__":
    main()