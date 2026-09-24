#!/usr/bin/env python3
"""
ir_error_stats.py

Computes MAE, RMSE, and max absolute deviation between calculated and
experimental IR band positions, from a CSV of consistently assigned
pairs (calculated cm-1, experimental cm-1) -- e.g., Table III of the
draft. Band matching/assignment is a human judgment call already made
in that table; this script only computes the error statistics, it
does not auto-match peaks.

Usage:
    python ~/scripts/ir_error_stats.py --pairs assigned_bands.csv
"""

import argparse
import csv
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Compute MAE/RMSE/max deviation for assigned IR bands.")
    parser.add_argument("--pairs", required=True,
                         help="CSV with columns: calculated_cm1,experimental_cm1,assignment")
    args = parser.parse_args()

    calc, exp, labels = [], [], []
    with open(args.pairs, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            calc.append(float(row["calculated_cm1"]))
            exp.append(float(row["experimental_cm1"]))
            labels.append(row.get("assignment", ""))

    calc = np.array(calc)
    exp = np.array(exp)
    diff = calc - exp

    mae = np.mean(np.abs(diff))
    rmse = np.sqrt(np.mean(diff ** 2))
    max_dev_idx = np.argmax(np.abs(diff))

    print(f"N bands: {len(calc)}")
    print(f"MAE:  {mae:.2f} cm-1")
    print(f"RMSE: {rmse:.2f} cm-1")
    print(f"Max absolute deviation: {np.abs(diff[max_dev_idx]):.2f} cm-1 "
          f"(calc {calc[max_dev_idx]:.1f}, exp {exp[max_dev_idx]:.1f}, "
          f"{labels[max_dev_idx]})")


if __name__ == "__main__":
    main()