#!/usr/bin/env python3
"""
plot_s1_scan.py

Relaxed S1 scan along the CNNC dihedral (output of parse_s1_scan.py):
  <output>_gap    : S1-S0 gap vs CNNC for one or more functionals
  <output>_energy : S0 and S1 energies (relative to the lowest S0) vs CNNC (first --scan only)
Only converged points (converged == 1) are drawn.

Usage:
    python ~/scripts/veusz/plot_s1_scan.py --scan CAM-B3LYP:s1_scan_camb3lyp.dat \
        --scan B3LYP:s1_scan_b3lyp.dat --output s1_scan
"""

import argparse

import pandas as pd

from _veusz_helpers import add_xy, finish, new_graph

HARTREE_TO_EV = 27.211386245988
COLORS = ["#d62728", "#1f77b4", "#2ca02c"]


def load(path):
    d = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None,
                    names=["target", "actual", "conv", "e_s0", "gap", "nm", "f", "orb", "w",
                           "ngap", "norb", "nw", "note"], na_values=["NA"])
    d = d[d.conv == 1].copy()
    for c in ("target", "e_s0", "gap", "f"):
        d[c] = pd.to_numeric(d[c])
    return d.sort_values("target")


def main():
    ap = argparse.ArgumentParser(description="Plot relaxed S1 scan (gap and energies).")
    ap.add_argument("--scan", action="append", required=True, help="LABEL:path")
    ap.add_argument("--output", default="s1_scan")
    args = ap.parse_args()

    scans = [(s.split(":", 1)[0], load(s.split(":", 1)[1])) for s in args.scan]

    g = new_graph("s1_gap", xlabel="CNNC dihedral (deg)", ylabel="S1 - S0 gap (eV)",
                  xmin=0, xmax=180, ymin=0)
    for i, (lab, d) in enumerate(scans):
        add_xy(g, f"gap{i}", d.target, d.gap, color=COLORS[i % 3], marker="circle", size="4pt",
               key=f"{lab} ({len(d)} converged pts)")
    g.Set("x/MajorTicks/number", 7)
    finish(g, f"{args.output}_gap", key_pos=("left", "bottom"))

    lab, d = scans[0]
    e0 = d.e_s0.min()
    s0 = (d.e_s0 - e0) * HARTREE_TO_EV
    s1 = s0 + d.gap
    g = new_graph("s1_energy", xlabel="CNNC dihedral (deg)",
                  ylabel="Energy relative to lowest S0 (eV)", xmin=0, xmax=180,
                  ymin=-0.2, ymax=float(s1.max()) + 0.9)
    add_xy(g, "s0", d.target, s0, color="black", marker="circle", size="4pt", key=f"S0 ({lab})")
    add_xy(g, "s1", d.target, s1, color="#d62728", marker="square", size="4pt",
           key=f"S1 ({lab}, S1 relaxed)")
    finish(g, f"{args.output}_energy", key_pos=("left", "top"))


if __name__ == "__main__":
    main()
