#!/usr/bin/env python3
"""
parse_bader.py

Parses the ACF.dat output from the Henkelman-group Bader code and
combines it with atomic symbols from a Turbomole coord file, computing
net Bader charges (all-electron reference: net_charge = Z - population)
and exporting a Veusz-importable dataset.

Usage:
    python ~/scripts/parse_bader.py --acf ACF.dat --coord coord \
        --output bader_dataset.dat
"""

import argparse
import numpy as np

BOHR_TO_ANGSTROM = 0.529177210903

# Atomic numbers (Z) for all-electron reference. Extend as needed.
ATOMIC_NUMBER = {
    "h": 1, "he": 2, "li": 3, "be": 4, "b": 5, "c": 6, "n": 7, "o": 8,
    "f": 9, "ne": 10, "na": 11, "mg": 12, "al": 13, "si": 14, "p": 15,
    "s": 16, "cl": 17, "ar": 18, "k": 19, "ca": 20,
}


def read_turbomole_coord_symbols(coord_path):
    """Read element symbols (1-indexed order) from a Turbomole $coord file."""
    symbols = []
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
                symbols.append(stripped.split()[3].lower())
    return symbols


def parse_acf(acf_path):
    """Parse ACF.dat: columns are # X Y Z CHARGE MIN DIST ATOMIC VOL."""
    charges = []
    with open(acf_path, "r") as f:
        lines = f.readlines()

    for line in lines:
        parts = line.split()
        if len(parts) >= 5 and parts[0].isdigit():
            charge = float(parts[4])
            charges.append(charge)

    return np.array(charges)


def main():
    parser = argparse.ArgumentParser(
        description="Parse Bader ACF.dat into a Veusz-importable dataset, "
                    "computing net charges for an all-electron (Turbomole) reference.")
    parser.add_argument("--acf", required=True, help="Path to ACF.dat")
    parser.add_argument("--coord", required=True, help="Path to Turbomole coord file")
    parser.add_argument("--output", default="bader_dataset.dat",
                         help="Output .dat filename (Veusz descriptor format)")
    args = parser.parse_args()

    symbols = read_turbomole_coord_symbols(args.coord)
    bader_populations = parse_acf(args.acf)

    if len(symbols) != len(bader_populations):
        print(f"WARNING: atom count mismatch - coord has {len(symbols)} atoms, "
              f"ACF.dat has {len(bader_populations)} entries. Check both files before proceeding.")

    unknown_elements = {s for s in symbols if s not in ATOMIC_NUMBER}
    if unknown_elements:
        print(f"WARNING: unknown element(s) not in ATOMIC_NUMBER table: {unknown_elements}. "
              f"Add their atomic numbers to the ATOMIC_NUMBER dict before trusting net charges.")

    with open(args.output, "w") as f:
        f.write("descriptor atom_index,element,bader_population,net_charge\n")
        for i, (symbol, population) in enumerate(zip(symbols, bader_populations), start=1):
            z = ATOMIC_NUMBER.get(symbol)
            if z is None:
                net_charge_str = "nan"
            else:
                net_charge = z - population
                net_charge_str = f"{net_charge:.6f}"
            f.write(f"{i} {symbol} {population:.6f} {net_charge_str}\n")

    print(f"Veusz dataset saved -> {args.output}")

    total_net_charge = sum(
        ATOMIC_NUMBER.get(s, 0) - p for s, p in zip(symbols, bader_populations)
        if s in ATOMIC_NUMBER
    )
    print(f"Sum of net charges (should be close to 0 for a neutral molecule): {total_net_charge:.4f}")


if __name__ == "__main__":
    main()