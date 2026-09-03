#!/usr/bin/env python3
"""
compare_charges.py

Cross-checks Bader net charges against Mulliken atomic charges,
atom by atom. Reports correlation, sign agreement, and largest
discrepancies -- intended to validate whether the Bader integration
error (~0.5 electrons, identified earlier) is distorting the chemical
interpretation, or just shifting the absolute scale.

Usage:
    python ~/scripts/compare_charges.py --bader bader_dataset.dat \
        --mulliken mulliken_dataset.dat --output charge_comparison.dat
"""

import argparse
import numpy as np


def read_descriptor_dat(path, columns):
    """
    Reads a Veusz-style .dat file with a 'descriptor' header line,
    extracting only the requested columns by name. Handles files that
    mix text columns (e.g. element symbols) with numeric ones by
    parsing line-by-line instead of np.loadtxt (which requires all
    columns to be numeric).
    """
    with open(path, "r") as f:
        header = f.readline()
        names = header.strip().replace("descriptor ", "").split(",")
        indices = [names.index(c) for c in columns]

        rows = []
        for line in f:
            parts = line.split()
            if len(parts) != len(names):
                continue
            rows.append([float(parts[i]) for i in indices])

    data = np.array(rows)
    return {c: data[:, i] for i, c in enumerate(columns)}


def main():
    parser = argparse.ArgumentParser(description="Cross-check Bader vs Mulliken atomic charges.")
    parser.add_argument("--bader", required=True, help="Path to bader_dataset.dat")
    parser.add_argument("--mulliken", required=True, help="Path to mulliken_dataset.dat")
    parser.add_argument("--output", default="charge_comparison.dat",
                         help="Output .dat filename (Veusz descriptor format)")
    args = parser.parse_args()

    bader = read_descriptor_dat(args.bader, ["atom_index", "net_charge"])
    mulliken = read_descriptor_dat(args.mulliken, ["atom_index", "mulliken_charge"])

    if not np.array_equal(bader["atom_index"], mulliken["atom_index"]):
        raise ValueError("Atom index columns don't match between the two files -- "
                          "check that both datasets cover the same 46 atoms in the same order.")

    bader_q = bader["net_charge"]
    mulliken_q = mulliken["mulliken_charge"]

    correlation = np.corrcoef(bader_q, mulliken_q)[0, 1]
    sign_agreement = np.mean(np.sign(bader_q) == np.sign(mulliken_q)) * 100
    diff = bader_q - mulliken_q

    print(f"Pearson correlation (Bader vs Mulliken): {correlation:.4f}")
    print(f"Sign agreement across all 46 atoms: {sign_agreement:.1f}%")
    print(f"Mean difference (Bader - Mulliken): {diff.mean():.4f}")
    print(f"Std dev of difference: {diff.std():.4f}")
    print()
    print("5 atoms with the LARGEST disagreement (Bader - Mulliken, absolute):")
    order = np.argsort(-np.abs(diff))
    for i in order[:5]:
        idx = int(bader["atom_index"][i])
        print(f"  Atom {idx}: Bader={bader_q[i]:.3f}, Mulliken={mulliken_q[i]:.3f}, "
              f"diff={diff[i]:.3f}")

    with open(args.output, "w") as f:
        f.write("descriptor atom_index,bader_net_charge,mulliken_charge\n")
        for i in range(len(bader_q)):
            f.write(f"{int(bader['atom_index'][i])} {bader_q[i]:.6f} {mulliken_q[i]:.6f}\n")

    print(f"\nCombined dataset saved -> {args.output}")


if __name__ == "__main__":
    main()