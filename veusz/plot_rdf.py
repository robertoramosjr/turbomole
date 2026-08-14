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
import numpy as np
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