#!/usr/bin/env python3
"""
plot_dos_fragments.py

Plots TDOS + fragment PDOS (named fragments, e.g. ring_1 / ring_2 /
azo_bridge) from parse_dos_fragments.py, with HOMO/LUMO dashed lines taken
from Multiwfn's orginfo.txt (or --homo/--lumo).

Usage:
    python ~/scripts/veusz/plot_dos_fragments.py --data dos_becke_fragments.dat \
        --orginfo orginfo.txt --output dos_becke_fragments_figure --emin -20 --emax 10
"""

import argparse

import numpy as np

from _veusz_helpers import add_vline, add_xy, finish, new_graph

HARTREE_TO_EV = 27.211386245988
COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]


def homo_lumo_from_orginfo(path):
    homo = lumo = None
    with open(path) as f:
        next(f)
        for line in f:
            p = line.split()
            if len(p) < 2:
                continue
            e, occ = float(p[0]), float(p[1])
            if occ > 0:
                homo = e
            elif lumo is None:
                lumo = e
    return homo * HARTREE_TO_EV, lumo * HARTREE_TO_EV


def main():
    ap = argparse.ArgumentParser(description="Plot TDOS + fragment PDOS in Veusz.")
    ap.add_argument("--data", required=True)
    ap.add_argument("--orginfo")
    ap.add_argument("--homo", type=float)
    ap.add_argument("--lumo", type=float)
    ap.add_argument("--emin", type=float, default=-20.0)
    ap.add_argument("--emax", type=float, default=10.0)
    ap.add_argument("--title", default="")
    ap.add_argument("--output", default="dos_fragments_figure")
    args = ap.parse_args()

    with open(args.data) as f:
        names = f.readline().split()[1].split(",")
    data = np.loadtxt(args.data, skiprows=1)
    energy = data[:, 0]
    win = (energy >= args.emin) & (energy <= args.emax)
    ymax = float(data[win, 1:].max()) * 1.5

    if args.homo is not None and args.lumo is not None:
        homo, lumo = args.homo, args.lumo
    elif args.orginfo:
        homo, lumo = homo_lumo_from_orginfo(args.orginfo)
    else:
        ap.error("give --orginfo or --homo/--lumo")

    g = new_graph("dos_fragments", xlabel="Energy (eV)", ylabel="DOS (states/eV, arb.)",
                  xmin=args.emin, xmax=args.emax, ymin=0, ymax=ymax)
    add_xy(g, "tdos", energy, data[:, 1], color="black", width="1.6pt", key=names[1])
    for i, nm in enumerate(names[2:]):
        add_xy(g, f"frag{i}", energy, data[:, 2 + i], color=COLORS[i % len(COLORS)],
               width="1.2pt", key=nm.replace("pDOS_", "PDOS "))
    add_vline(g, "homo", homo, 0, ymax, key=f"HOMO {homo:.2f} eV")
    add_vline(g, "lumo", lumo, 0, ymax, color="darkgray", style="dotted",
              key=f"LUMO {lumo:.2f} eV (gap {lumo - homo:.2f} eV)")
    finish(g, args.output, key_pos=("right", "top"))


if __name__ == "__main__":
    main()
