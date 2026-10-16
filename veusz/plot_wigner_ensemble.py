#!/usr/bin/env python3
"""
plot_wigner_ensemble.py

From a parse_wigner_ic_dump.py CSV (1000 Wigner points x N states):
  <output>_hist  : histogram of the S1 and S2 vertical energies, with the laser
                   windows (center = laser energy + functional bias, +- halfwidth)
  <output>_scatter: oscillator strength vs energy (log y) for S1 and S2

Usage:
    python ~/scripts/veusz/plot_wigner_ensemble.py \
        --dump wigner_1000_states_dump.csv --output wigner_ensemble \
        --laser 325:S2:0.147 --laser 425:S1:-0.109 --halfwidth 0.10
"""

import argparse

import numpy as np
import pandas as pd

from _veusz_helpers import add_vline, add_xy, finish, new_graph

HC = 1239.84198  # eV nm


def main():
    ap = argparse.ArgumentParser(description="Plot Wigner-ensemble excitation energies.")
    ap.add_argument("--dump", required=True)
    ap.add_argument("--dump2", help="Optional second dump (e.g. B3LYP) overlaid in the histogram")
    ap.add_argument("--label", default="CAM-B3LYP")
    ap.add_argument("--label2", default="B3LYP")
    ap.add_argument("--laser", action="append", default=[],
                    help="lambda_nm:S1|S2:bias_eV (bias = functional error vs experiment; "
                         "model center = laser_eV + bias)")
    ap.add_argument("--halfwidth", type=float, default=0.10)
    ap.add_argument("--bin", type=float, default=0.05)
    ap.add_argument("--output", default="wigner_ensemble")
    args = ap.parse_args()

    d = pd.read_csv(args.dump)
    d2 = pd.read_csv(args.dump2) if args.dump2 else None
    edges = np.arange(1.5, 6.0 + args.bin, args.bin)
    mid = 0.5 * (edges[1:] + edges[:-1])

    def step(h):
        xs, ys = [], []
        for lo, hi, v in zip(edges[:-1], edges[1:], h):
            xs += [lo, hi]
            ys += [v, v]
        return xs, ys

    g = new_graph("wigner_hist", xlabel="Vertical excitation energy (eV)",
                  ylabel="Wigner points per bin", xmin=1.8, xmax=5.2, ymin=0)
    hmax = 0
    for s, col in ((1, "#1f77b4"), (2, "#d62728")):
        h, _ = np.histogram(d[d.state == s].energy_eV, bins=edges)
        hmax = max(hmax, h.max())
        xs, ys = step(h)
        add_xy(g, f"s{s}", xs, ys, color=col, width="1.4pt", key=f"S{s} {args.label}")
        if d2 is not None:
            h2, _ = np.histogram(d2[d2.state == s].energy_eV, bins=edges)
            xs, ys = step(h2)
            add_xy(g, f"s{s}b", xs, ys, color=col, width="1pt", style="dashed",
                   key=f"S{s} {args.label2}")
    g.Set("y/max", float(hmax * 1.75))
    for spec in args.laser:
        nm, st, bias = spec.split(":")
        c = HC / float(nm) + float(bias)
        add_vline(g, f"laser{nm}", c, 0, hmax * 1.75, color="black", style="dashed",
                  key=f"{nm} nm -> {c:.3f} eV ({st}, bias {float(bias):+.3f})")
        for sgn in (-1, 1):
            add_vline(g, f"win{nm}_{sgn}", c + sgn * args.halfwidth, 0, hmax * 1.75,
                      color="gray", style="dotted")
    finish(g, f"{args.output}_hist", key_pos=("right", "top"))

    g = new_graph("wigner_scatter", xlabel="Vertical excitation energy (eV)",
                  ylabel="Oscillator strength (length)", xmin=1.8, xmax=5.2, ymin=1e-6, ymax=2)
    g.Set("y/log", True)
    for s, col in ((1, "#1f77b4"), (2, "#d62728")):
        sub = d[(d.state == s) & (d.osc_len > 0)]
        add_xy(g, f"p{s}", sub.energy_eV, sub.osc_len, color=col, marker="circle", size="1.5pt",
               line=False, key=f"S{s}")
    for spec in args.laser:
        nm, st, bias = spec.split(":")
        c = HC / float(nm) + float(bias)
        add_vline(g, f"laser{nm}", c, 1e-6, 2, color="black", style="dashed",
                  key=f"{nm} nm ({c:.3f} eV)")
    finish(g, f"{args.output}_scatter", key_pos=("left", "top"))


if __name__ == "__main__":
    main()
