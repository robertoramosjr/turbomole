#!/usr/bin/env python3
"""
plot_svp_tda_correlation.py

Shift of the S1/S2 vertical energies between def2-SVP/TDA and def2-TZVP/RPA
for the same Wigner points, plotted against the TZVP/RPA energy. A flat
cloud = rigid offset (S1); a sloped cloud = non-rigid (S2).

Usage:
    python ~/scripts/veusz/plot_svp_tda_correlation.py \
        --csv svp_tda_shift_comparison.csv --output svp_tda_correlation
"""

import argparse

import numpy as np
import pandas as pd

from _veusz_helpers import add_xy, finish, new_graph


def main():
    ap = argparse.ArgumentParser(description="Plot SVP/TDA - TZVP/RPA shifts.")
    ap.add_argument("--csv", required=True)
    ap.add_argument("--output", default="svp_tda_correlation")
    args = ap.parse_args()

    d = pd.read_csv(args.csv)
    g = new_graph("svp_tda", xlabel="E(TZVP/RPA) (eV)",
                  ylabel="E(SVP/TDA) - E(TZVP/RPA) (eV)", xmin=1.8, xmax=5.0)
    for s, col, lab in ((1, "#1f77b4", "S1 (n->pi*)"), (2, "#d62728", "S2 (pi->pi*)")):
        x = d[f"s{s}_tzvp_rpa_eV"].to_numpy()
        y = d[f"s{s}_shift_eV"].astype(float).to_numpy()
        add_xy(g, f"s{s}", x, y, color=col, marker="circle", size="4pt", line=False,
               key=f"{lab}: {y.mean():+.3f} +- {y.std(ddof=1):.3f} eV")
        m, b = np.polyfit(x, y, 1)
        xs = np.array([x.min(), x.max()])
        add_xy(g, f"fit{s}", xs, m * xs + b, color=col, width="1pt", style="dashed",
               key=f"fit S{s}: r = {np.corrcoef(x, y)[0, 1]:+.2f}")
    add_xy(g, "zero", [1.8, 5.0], [0, 0], color="gray", width="0.8pt", style="dotted")
    finish(g, args.output, key_pos=("left", "top"))


if __name__ == "__main__":
    main()
