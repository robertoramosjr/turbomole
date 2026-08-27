#!/usr/bin/env python3
"""
plot_rdf.py

Builds an N-panel grid of RDF plots from a Veusz-style .dat file
(descriptor header), one panel per --panel argument. Fully
parameterized via argparse -- no need to edit the script per
molecule; just pass the right --panel flags on the command line.

Usage:
    python ~/scripts/veusz/plot_rdf.py --data rdf_dataset.dat \
        --output rdf_figure \
        --panel gr_carboxylic_acid "Carboxylic Acid" black a \
        --panel gr_phenolic_ring "Phenolic Ring" red b \
        --panel gr_prenyl_1 "Prenyl 1" green c \
        --panel gr_prenyl_2 "Prenyl 2" blue d
"""

import argparse
import os

import numpy as np

if not os.environ.get("DISPLAY"):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import veusz.embed as veusz


def compute_ymax(values, percentile=95, headroom=1.3):
    nonzero = values[values > 0]
    if nonzero.size == 0:
        return 1.0
    return float(np.percentile(nonzero, percentile) * headroom)


def add_rdf_panel(g, name, dataset_name, title, color, letter, ymax):
    g.Add("graph", name=name)
    g.To(name)

    g.Set("x/label", "r (Angstrom)")
    g.Set("x/min", 0)
    g.Set("x/max", 14)

    g.Set("y/label", "g(r)")
    g.Set("y/min", 0)
    g.Set("y/max", ymax)

    xy_name = f"xy_{dataset_name}"
    g.Add("xy", name=xy_name)
    g.Set(f"{xy_name}/xData", "r_angstrom")
    g.Set(f"{xy_name}/yData", dataset_name)
    g.Set(f"{xy_name}/marker", "none")
    g.Set(f"{xy_name}/PlotLine/width", "1pt")
    g.Set(f"{xy_name}/PlotLine/color", color)

    g.Add("label", name="title_label", label=title)
    g.Set("title_label/xPos", 0.5)
    g.Set("title_label/yPos", 1.05)
    g.Set("title_label/alignHorz", "centre")

    g.Add("label", name="letter_label", label=f"{letter})")

    g.To("..")


def read_header_names(path):
    with open(path, "r") as f:
        first = f.readline().strip()
    return first[len("descriptor"):].strip().split(",")


def main():
    parser = argparse.ArgumentParser(description="Plot an N-panel grid of RDFs in Veusz.")
    parser.add_argument("--data", required=True, help="Path to rdf_dataset.dat")
    parser.add_argument("--output", default="rdf_figure", help="Output filename prefix")
    parser.add_argument(
        "--panel", action="append", nargs=4, required=True,
        metavar=("DATASET", "TITLE", "COLOR", "LETTER"),
        help="Add a panel: dataset column name, title, line color, panel letter "
             "(may be repeated, one --panel per RDF curve)",
    )
    args = parser.parse_args()

    names = read_header_names(args.data)
    raw = np.loadtxt(args.data, skiprows=1)

    g = veusz.Embedded("rdf_figure", hidden=True)

    g.To(g.Add("page", name="page1"))
    g.To(g.Add("grid", name="grid1", columns=2))

    g.ImportFile(args.data, "", linked=True)

    for dataset_name, title, color, letter in args.panel:
        col = names.index(dataset_name)
        ymax = compute_ymax(raw[:, col])
        add_rdf_panel(g, f"graph_{dataset_name}", dataset_name, title, color, letter, ymax)

    g.Save(f"{args.output}.vsz")
    g.Export(f"{args.output}.pdf")
    print(f"Saved -> {args.output}.vsz and {args.output}.pdf")
    g.Close()


if __name__ == "__main__":
    main()