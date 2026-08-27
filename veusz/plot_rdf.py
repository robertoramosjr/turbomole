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
    g.Set("letter_label/xPos", 0.05)
    g.Set("letter_label/yPos", 0.92)
    g.Set("letter_label/Text/size", "14pt")
    g.Set("letter_label/Text/bold", True)

    g.To("..")


def load_columns(path):
    with open(path) as f:
        header = f.readline().strip()
    names = header.split(None, 1)[1].split(",")
    data = np.loadtxt(path, skiprows=1)
    return {name: data[:, i] for i, name in enumerate(names)}


def main():
    parser = argparse.ArgumentParser(description="Plot an N-panel RDF grid in Veusz.")
    parser.add_argument("--data", required=True, help="Path to rdf_dataset.dat")
    parser.add_argument("--output", default="rdf_figure", help="Output filename prefix")
    parser.add_argument("--panel", nargs=4, action="append", required=True,
                         metavar=("DATASET", "TITLE", "COLOR", "LETTER"),
                         help="Panel spec: dataset_name title color letter. Repeat once per panel.")
    parser.add_argument("--columns", type=int, default=2, help="Number of grid columns")
    args = parser.parse_args()

    columns_data = load_columns(args.data)

    g = veusz.Embedded("rdf_figure", hidden=True)

    g.To(g.Add("page", name="page1"))
    g.Set("width", "16cm")
    g.Set("height", "13cm")
    g.To(g.Add("grid", name="grid1", columns=args.columns))

    g.ImportFile(args.data, "", linked=True)

    for dataset_name, title, color, letter in args.panel:
        ymax = compute_ymax(columns_data[dataset_name])
        add_rdf_panel(g, f"graph_{dataset_name}", dataset_name, title, color, letter, ymax)

    g.Save(f"{args.output}.vsz")
    g.Export(f"{args.output}.pdf")
    print(f"Saved -> {args.output}.vsz and {args.output}.pdf")
    g.Close()


if __name__ == "__main__":
    main()