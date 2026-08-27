#!/usr/bin/env python3
"""
plot_becke_pdos.py

Plots the Becke-partitioned TDOS and per-element PDOS (O/C/H) produced by
parse_becke_pdos.py, all overlaid in a single graph, with the HOMO and
LUMO energies marked as vertical dashed lines and the gap between them
annotated on the plot.

HOMO/LUMO are not re-derived from the (broadened) DOS curve -- the curve
doesn't dip back to zero at this broadening width, so there's no reliable
zero-crossing to detect. Instead they're read from Multiwfn's own
orginfo.txt (written alongside DOS_curve.txt: one row per orbital, energy
in a.u. + occupation), by default -- HOMO is the highest-energy orbital
with occupation > 0, LUMO the lowest-energy orbital with occupation == 0.
Cross-checked against Turbomole's own report for this dataset (proper.out,
"HOMO-LUMO Separation"): agrees to better than 0.001 eV. Pass --homo/--lumo
directly instead if you already have them and don't want to rely on
orginfo.txt.

Usage:
    python ~/scripts/veusz/plot_becke_pdos.py --data dos_becke_dataset.dat \
        --orginfo orginfo.txt --output becke_pdos_figure

    # or, with HOMO/LUMO given directly (e.g. from
    # grep -A 4 "HOMO-LUMO Separation" proper.out):
    python ~/scripts/veusz/plot_becke_pdos.py --data dos_becke_dataset.dat \
        --homo -5.79375 --lumo -2.10956 --output becke_pdos_figure
"""

import argparse
import os

import numpy as np

if not os.environ.get("DISPLAY"):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import veusz.embed as veusz

HARTREE_TO_EV = 27.211386245988


def homo_lumo_from_orginfo(path):
    """Derive HOMO/LUMO (eV) from Multiwfn's orginfo.txt: one row per
    orbital in ascending energy order, columns are energy (a.u.),
    occupation, weight, broadening width. HOMO is the highest-energy
    orbital with occupation > 0, LUMO the lowest-energy one with
    occupation == 0."""
    homo_au = None
    lumo_au = None
    with open(path, "r") as f:
        next(f)  # header: n_orbitals, n_spins
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            energy_au, occ = float(parts[0]), float(parts[1])
            if occ > 0:
                homo_au = energy_au
            elif lumo_au is None:
                lumo_au = energy_au
    if homo_au is None or lumo_au is None:
        raise ValueError(
            f"Could not find both an occupied and an unoccupied orbital "
            f"in {path} -- check it's a Multiwfn orginfo.txt with an "
            f"occupation column and orbitals in ascending energy order.")
    return homo_au * HARTREE_TO_EV, lumo_au * HARTREE_TO_EV


def add_vertical_line(graph, name, x_ev, color):
    line = graph.Add("line", name=name)
    line.positioning.val = "axes"
    line.mode.val = "point-to-point"
    line.xPos.val = [x_ev]
    line.yPos.val = [0]
    line.xPos2.val = [x_ev]
    line.yPos2.val = [1e6]  # clipped to the graph's actual y range
    line.Line.color.val = color
    line.Line.style.val = "dashed"
    line.Line.width.val = "1pt"
    return line


def main():
    parser = argparse.ArgumentParser(
        description="Plot Becke-partitioned TDOS + per-element PDOS in a single "
                    "graph, with HOMO/LUMO markers and the gap annotated.")
    parser.add_argument("--data", required=True, help="Path to dos_becke_dataset.dat")
    parser.add_argument("--output", default="becke_pdos_figure", help="Output filename prefix")
    parser.add_argument("--emin", type=float, default=-20.0)
    parser.add_argument("--emax", type=float, default=10.0)
    parser.add_argument("--orginfo",
                         help="Path to Multiwfn's orginfo.txt, used to derive HOMO/LUMO "
                              "automatically (default, unless --homo/--lumo are given)")
    parser.add_argument("--homo", type=float,
                         help="HOMO energy in eV -- overrides --orginfo if given "
                              "(e.g. from 'grep -A 4 \"HOMO-LUMO Separation\" proper.out')")
    parser.add_argument("--lumo", type=float,
                         help="LUMO energy in eV -- overrides --orginfo if given")
    args = parser.parse_args()

    if args.homo is not None and args.lumo is not None:
        homo, lumo = args.homo, args.lumo
    elif args.orginfo:
        homo, lumo = homo_lumo_from_orginfo(args.orginfo)
        print(f"HOMO/LUMO derived from {args.orginfo}")
    else:
        parser.error("Provide either --homo and --lumo directly, or --orginfo "
                      "to derive them automatically from Multiwfn's output.")

    gap = lumo - homo

    energy, tdos, pdos_o, pdos_c, pdos_h = np.loadtxt(
        args.data, skiprows=1, unpack=True)
    in_window = (energy >= args.emin) & (energy <= args.emax)
    ymax = float(max(tdos[in_window].max(), pdos_o[in_window].max(),
                      pdos_c[in_window].max(), pdos_h[in_window].max()))

    g = veusz.Embedded("becke_pdos_figure", hidden=True)
    g.EnableToolbar(False)

    page = g.Root.Add("page", name="page1")
    graph = page.Add("graph", name="graph1")

    g.To(graph.path)
    g.ImportFile(args.data, "", linked=True)

    # Set properties on the graph's default x/y axes rather than adding new
    # ones -- graph.Add("axis", ...) creates a duplicate ghost axis on top
    # of the auto-created default, doubling up the tick labels.
    graph.x.label.val = "Energy (eV)"
    graph.x.min.val = args.emin
    graph.x.max.val = args.emax

    graph.y.label.val = "DOS (arb. units)"
    graph.y.min.val = 0
    # Explicit headroom above the tallest curve so the HOMO/LUMO/gap labels
    # (placed in axis coordinates) aren't clipped by an auto-ranged top.
    graph.y.max.val = ymax * 1.2

    curves = [
        ("xy_tdos", "energy_eV", "TDOS", "Total", "black", "2pt"),
        ("xy_pdos_o", "energy_eV", "pDOS_O", "O", "red", "1pt"),
        ("xy_pdos_c", "energy_eV", "pDOS_C", "C", "blue", "1pt"),
        ("xy_pdos_h", "energy_eV", "pDOS_H", "H", "green", "1pt"),
    ]
    for name, xdata, ydata, key, color, width in curves:
        xy = graph.Add("xy", name=name)
        xy.xData.val = xdata
        xy.yData.val = ydata
        xy.marker.val = "none"
        xy.PlotLine.color.val = color
        xy.PlotLine.width.val = width
        xy.key.val = key

    add_vertical_line(graph, "homo_line", homo, "grey")
    add_vertical_line(graph, "lumo_line", lumo, "grey")

    homo_label = graph.Add("label", name="homo_label")
    homo_label.positioning.val = "axes"
    homo_label.xPos.val = [homo]
    homo_label.yPos.val = [ymax * 1.02]
    homo_label.label.val = "HOMO"
    homo_label.alignHorz.val = "right"

    lumo_label = graph.Add("label", name="lumo_label")
    lumo_label.positioning.val = "axes"
    lumo_label.xPos.val = [lumo]
    lumo_label.yPos.val = [ymax * 1.02]
    lumo_label.label.val = "LUMO"
    lumo_label.alignHorz.val = "left"

    gap_label = graph.Add("label", name="gap_label")
    gap_label.positioning.val = "axes"
    gap_label.xPos.val = [(homo + lumo) / 2]
    gap_label.yPos.val = [ymax * 1.12]
    gap_label.label.val = f"Gap = {gap:.2f} eV"
    gap_label.alignHorz.val = "centre"

    graph.Add("key", name="key1")
    g.Set("key1/Text/size", "9pt")

    g.Save(f"{args.output}.vsz")
    g.Export(f"{args.output}.pdf")
    print(f"HOMO = {homo:.4f} eV, LUMO = {lumo:.4f} eV, gap = {gap:.4f} eV")
    print(f"Saved -> {args.output}.vsz and {args.output}.pdf")
    g.Close()


if __name__ == "__main__":
    main()
