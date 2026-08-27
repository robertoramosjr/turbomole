#!/usr/bin/env python3
"""
plot_dos.py

Builds a 2x2 grid: (a) TDOS with s/p/d contributions, (b) JDOS-style
profile (Gaussian-broadened TD-DFT oscillator strengths), (c) pDOS-C/H/O,
(d) reserved. Uses the classic (name-based) Veusz embedding API, with
default graph axes only (no Add('axis', ...) calls) to avoid the
ghost-axis duplication seen in earlier versions.

Usage:
    python ~/scripts/veusz/plot_dos.py --data dos_dataset.dat \
        --jdos jdos_dataset.dat --output dos_figure --emin -20 --emax 10
"""

import argparse
import os

if not os.environ.get("DISPLAY"):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import veusz.embed as veusz


def add_multiline_graph(g, name, title, curves, xlabel, ylabel, emin, emax, xdata="energy_eV"):
    g.Add("graph", name=name)
    g.To(name)

    g.Set("x/label", xlabel)
    g.Set("x/min", emin)
    g.Set("x/max", emax)

    g.Set("y/label", ylabel)
    g.Set("y/autoRange", "next")

    for dataset_name, legend_label, color in curves:
        xy_name = f"xy_{dataset_name}"
        g.Add("xy", name=xy_name)
        g.Set(f"{xy_name}/xData", xdata)
        g.Set(f"{xy_name}/yData", dataset_name)
        g.Set(f"{xy_name}/marker", "none")
        g.Set(f"{xy_name}/PlotLine/color", color)
        g.Set(f"{xy_name}/key", legend_label)

    g.Add("label", name="title_label", label=title)
    g.Set("title_label/xPos", 0.5)
    g.Set("title_label/yPos", 1.05)
    g.Set("title_label/alignHorz", "centre")

    g.Add("key", name="key1")
    g.Set("key1/Text/size", "8pt")

    g.To("..")


def add_placeholder_graph(g, name, text):
    g.Add("graph", name=name)
    g.To(name)
    g.Set("x/hide", True)
    g.Set("y/hide", True)
    g.Add("label", name="placeholder_label", label=text)
    g.Set("placeholder_label/xPos", 0.5)
    g.Set("placeholder_label/yPos", 0.5)
    g.Set("placeholder_label/alignHorz", "centre")
    g.To("..")


def main():
    parser = argparse.ArgumentParser(description="Plot TDOS/JDOS/PDOS panels in Veusz.")
    parser.add_argument("--data", required=True, help="Path to dos_dataset.dat")
    parser.add_argument("--jdos", required=True, help="Path to jdos_dataset.dat")
    parser.add_argument("--output", default="dos_figure", help="Output filename prefix")
    parser.add_argument("--emin", type=float, default=-20.0)
    parser.add_argument("--emax", type=float, default=10.0)
    args = parser.parse_args()

    g = veusz.Embedded("dos_figure", hidden=True)

    g.To(g.Add("page", name="page1"))
    g.To(g.Add("grid", name="grid1", columns=2))

    g.ImportFile(args.data, "", linked=True)
    g.ImportFile(args.jdos, "", linked=True)

    add_multiline_graph(
        g, "graph_tdos", "TDOS",
        [("TDOS", "Total", "black"), ("DOS_s", "s", "red"), ("DOS_p", "p", "blue")],
        "Energy (eV)", "DOS (states/eV)", args.emin, args.emax
    )

    # JDOS-style panel: broadened TD-DFT oscillator strengths, plotted over
    # its own native energy range (excited-state energies), not the DOS emin/emax.
    add_multiline_graph(
        g, "graph_jdos", "JDOS (TD-DFT, broadened)",
        [("jdos_intensity", "Broadened osc. strength", "black")],
        "Energy (eV)", "Intensity (arb. units)", 0, 10,
        xdata="jdos_energy_eV"
    )

    add_multiline_graph(
        g, "graph_pdos", "pDOS by element",
        [("pDOS_C", "C", "black"), ("pDOS_H", "H", "grey"), ("pDOS_O", "O", "red")],
        "Energy (eV)", "PDOS (states/eV)", args.emin, args.emax
    )

    add_placeholder_graph(g, "graph_reserved", "(reserved)")

    g.Save(f"{args.output}.vsz")
    g.Export(f"{args.output}.pdf")
    print(f"Saved -> {args.output}.vsz and {args.output}.pdf")
    g.Close()


if __name__ == "__main__":
    main()