#!/usr/bin/env python3
"""
group_average_charges.py

Computes group-averaged Bader net charges (mean and range) from
groups.json (functional group -> atom indices) and bader_dataset.dat
(per-atom net charges), restoring the group-level narrative
(carboxylic vs. phenolic vs. prenyl) for the paper text.

Usage:
    python ~/scripts/group_average_charges.py --groups groups.json \
        --bader bader_dataset.dat
"""

import argparse
import json
import numpy as np


def read_bader_charges(path):
    charges = {}
    with open(path, "r") as f:
        header = f.readline()
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            atom_index = int(parts[0])
            net_charge = float(parts[3])  # atom_index, element, bader_population, net_charge
            charges[atom_index] = net_charge
    return charges

def main():
    parser = argparse.ArgumentParser(description="Compute group-averaged Bader charges.")
    parser.add_argument("--groups", required=True, help="Path to groups.json")
    parser.add_argument("--bader", required=True, help="Path to bader_dataset.dat")
    args = parser.parse_args()

    with open(args.groups, "r") as f:
        groups = json.load(f)

    charges = read_bader_charges(args.bader)

    print(f"{'Group':<20} {'N':<6} {'Mean':<12} {'Mean|q|':<12} {'Min':<10} {'Max':<10} {'Range':<10}")
    for group_name, indices in groups.items():
        values = np.array([charges[i] for i in indices if i in charges])
        mean_abs = np.mean(np.abs(values))
        value_range = values.max() - values.min()
        print(f"{group_name:<20} {len(values):<6} {values.mean():<12.4f} {mean_abs:<12.4f} "
              f"{values.min():<10.4f} {values.max():<10.4f} {value_range:<10.4f}")


if __name__ == "__main__":
    main()