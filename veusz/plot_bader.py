#!/usr/bin/env python3
"""
plot_bader.py

Bar plot of net Bader charges per atom index, from bader_dataset.dat.
Mirrors Fig. 2 of Januario & Cabral (2026), as a bar chart instead of
charges annotated on the 3D structure (structure annotation is better
done in a molecular viewer, not Veusz).

Usage:
    python ~/scripts/plot_bader.py --data bader_dataset.dat --output bader_figure
"""

import argparse
import os

if not os.environ.get("DISPLAY"):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import veusz.embed as veusz


def main():
    parser = argparse.ArgumentParser(description="Plot net Bader charges in Veusz.")
    parser.add_argument("--data", required=True, help="Path to bader_dataset.dat")
    parser.add_argument("--output", default="bader_figure", help="Output filename prefix")
    args = parser.parse_args()

    g = veusz.Embedded("bader_figure", hidden=True)
    g.EnableToolbar(False)

    page = g.Root.Add("page", name="page1")
    graph = page.Add("graph", name="graph1")

    g.To(graph.path)
    g.ImportFile(args.data, "", linked=True)

    # Set properties on the graph's default x/y axes rather than adding new
    # ones -- graph.Add("axis", ...) creates a duplicate ghost axis on top
    # of the auto-created default, doubling up the tick labels.
    graph.x.label.val = "Atom index"
    graph.y.label.val = "Net Bader charge (e)"

    bar = graph.Add("bar", name="bar1")
    bar.lengths.val = ["net_charge"]
    bar.posn.val = "atom_index"

    g.Save(f"{args.output}.vsz")
    g.Export(f"{args.output}.pdf")
    print(f"Saved -> {args.output}.vsz and {args.output}.pdf")
    g.Close()


if __name__ == "__main__":
    main()