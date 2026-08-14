#!/usr/bin/env python3
"""
plot_bond_order.py

Bar plot of the top-N strongest atom-pair shared electron numbers
(SEN), a Mulliken-derived bond-strength descriptor used here as the
closest native Turbomole analogue to VASP+LOBSTER's COHP/ICOHP.

Usage:
    python ~/scripts/veusz/plot_bond_order.py --data bond_order_dataset.dat \
        --labels bond_order_dataset_labels.csv --top 20 --output bond_order_figure
"""

import argparse
import csv
import veusz.embed as veusz


def read_labels_csv(path, top_n):
    rows = []
    with open(path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "pair_index": int(row["pair_index"]),
                "pair_label": row["pair_label"],
                "sen": float(row["shared_electron_number"]),
            })
    rows.sort(key=lambda r: r["sen"], reverse=True)
    return rows[:top_n]


def main():
    parser = argparse.ArgumentParser(description="Plot top-N SEN bond-order pairs in Veusz.")
    parser.add_argument("--data", required=True, help="Path to bond_order_dataset.dat (unused directly; kept for consistency)")
    parser.add_argument("--labels", required=True, help="Path to bond_order_dataset_labels.csv")
    parser.add_argument("--top", type=int, default=20, help="Number of strongest pairs to display")
    parser.add_argument("--output", default="bond_order_figure", help="Output filename prefix")
    args = parser.parse_args()

    top_pairs = read_labels_csv(args.labels, args.top)

    print(f"Plotting top {len(top_pairs)} pairs by SEN:")
    for i, p in enumerate(top_pairs, start=1):
        print(f"  {i}. {p['pair_label']}: {p['sen']:.4f}")

    # Write a small ranked dataset (rank -> SEN) for Veusz to import.
    ranked_path = f"{args.output}_ranked.dat"
    with open(ranked_path, "w") as f:
        f.write("descriptor rank,sen_value\n")
        for i, p in enumerate(top_pairs, start=1):
            f.write(f"{i} {p['sen']:.6f}\n")

    g = veusz.Embedded("bond_order_figure", hidden=True)

    g.To(g.Add("page", name="page1"))
    g.To(g.Add("graph", name="graph1"))

    g.ImportFile(ranked_path, "", linked=True)

    g.Set("x/label", "Bond pair rank (see console output / labels CSV for atom identities)")
    g.Set("x/min", 0)
    g.Set("x/max", len(top_pairs) + 1)

    g.Set("y/label", "Shared electron number (SEN)")
    g.Set("y/min", 0)

    g.Add("bar", name="bar1")
    g.Set("bar1/lengths", ["sen_value"])
    g.Set("bar1/posn", "rank")

    g.Save(f"{args.output}.vsz")
    g.Export(f"{args.output}.pdf")
    print(f"Saved -> {args.output}.vsz and {args.output}.pdf")
    g.Close()


if __name__ == "__main__":
    main()