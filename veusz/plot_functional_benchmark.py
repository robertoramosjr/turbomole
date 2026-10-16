#!/usr/bin/env python3
"""
plot_functional_benchmark.py

Vertical-excitation benchmark (trans-azobenzene): signed error vs experiment
(n->pi*, pi->pi*) per functional, from functional_benchmark_trans.csv.

Usage:
    python ~/scripts/veusz/plot_functional_benchmark.py \
        --csv functional_benchmark_trans.csv --output functional_benchmark_figure
"""

import argparse
import csv

from _veusz_helpers import add_xy, finish, new_graph

ORDER = ["pbe0", "b3lyp", "hse06", "camb3lyp"]
LABEL = {"pbe0": "PBE0", "b3lyp": "B3LYP", "hse06": "HSE06", "camb3lyp": "CAM-B3LYP"}


def main():
    ap = argparse.ArgumentParser(description="Plot functional benchmark errors.")
    ap.add_argument("--csv", required=True)
    ap.add_argument("--output", default="functional_benchmark_figure")
    args = ap.parse_args()

    err = {"n->pi*": {}, "pi->pi*": {}}
    with open(args.csv) as f:
        for r in csv.DictReader(f):
            if r["reference_label"] in err:
                err[r["reference_label"]][r["run"]] = float(r["reference_error_eV"])
    runs = [r for r in ORDER if r in err["n->pi*"]]

    g = new_graph("functional_benchmark", xlabel="Functional (1 PBE0, 2 B3LYP, 3 HSE06, 4 CAM-B3LYP)",
                  ylabel="Error vs experiment (eV)", xmin=0.4, xmax=len(runs) + 0.6,
                  ymin=-0.4, ymax=0.3)
    xs = list(range(1, len(runs) + 1))
    add_xy(g, "zero", [0.4, len(runs) + 0.6], [0, 0], color="gray", width="0.8pt", style="dotted")
    add_xy(g, "npi", [x - 0.08 for x in xs], [err["n->pi*"][r] for r in runs], color="#1f77b4",
           marker="square", size="7pt", line=False, key="n->pi* (exp. 2.82 eV)")
    add_xy(g, "pipi", [x + 0.08 for x in xs], [err["pi->pi*"][r] for r in runs], color="#d62728",
           marker="circle", size="7pt", line=False, key="pi->pi* (exp. 3.90 eV)")
    # tick labels: the x axis is numeric; label each position with a text widget
    g.Set("x/MajorTicks/manualTicks", [float(x) for x in xs])
    for x, r in zip(xs, runs):
        g.Add("label", name=f"lab_{r}")
        g.Set(f"lab_{r}/label", LABEL[r])
        g.Set(f"lab_{r}/xPos", [(x - 0.4) / (len(runs) + 0.2)])
        g.Set(f"lab_{r}/yPos", [0.04])
        g.Set(f"lab_{r}/Text/size", "9pt")
    finish(g, args.output, key_pos=("left", "top"))


if __name__ == "__main__":
    main()
