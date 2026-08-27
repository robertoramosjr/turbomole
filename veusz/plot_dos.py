#!/usr/bin/env python3
"""
plot_dos.py

Builds a 2x2 grid: (a) TDOS with s/p/d contributions, (b) JDOS-style
profile (Gaussian-broadened TD-DFT oscillator strengths), (c) pDOS-C/H/<4th
element>, (d) reserved. Uses the classic (name-based) Veusz embedding API,
with default graph axes only (no Add('axis', ...) calls) to avoid the
ghost-axis duplication seen in earlier versions.

Optionally overlays the HOMO/LUMO levels (from parse_homo_lumo.py's
homo_lumo_dataset.dat) as vertical dashed lines on the TDOS panel, with
the gap value labelled -- otherwise that number only ever lived in a
.out file's terminal output.

Usage:
    python ~/scripts/veusz/plot_dos.py --data dos_dataset.dat \
        --jdos jdos_dataset.dat --output dos_figure --emin -20 --emax 10 \
        --element4-label O --homo-lumo homo_lumo_dataset.dat

    # azobenzeno (N no lugar de O):
    python ~/scripts/veusz/plot_dos.py --data dos_dataset.dat \
        --jdos jdos_dataset.dat --output dos_figure --emin -20 --emax 10 \
        --element4-label N --homo-lumo homo_lumo_dataset.dat
"""

import argparse
import os

if not os.environ.get("DISPLAY"):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import veusz.embed as veusz

# Drawn past any realistic DOS y-range so the vertical HOMO/LUMO lines
# fill the panel; Veusz clips xy data to the visible axis rectangle.
VLINE_Y_EXTENT = 1.0e4


def read_homo_lumo(path):
    with open(path, "r") as f:
        lines = [line for line in f if not line.startswith("descriptor")]
    homo_ev, lumo_ev, gap_ev = (float(x) for x in lines[0].split())
    return {"homo_eV": homo_ev, "lumo_eV": lumo_ev, "gap_eV": gap_ev}


def add_homo_lumo_lines(g, homo_lumo):
    g.SetData("homo_line_x", [homo_lumo["homo_eV"], homo_lumo["homo_eV"]])
    g.SetData("homo_line_y", [-VLINE_Y_EXTENT, VLINE_Y_EXTENT])
    g.SetData("lumo_line_x", [homo_lumo["lumo_eV"], homo_lumo["lumo_eV"]])
    g.SetData("lumo_line_y", [-VLINE_Y_EXTENT, VLINE_Y_EXTENT])

    for name, key in [("homo_line", "HOMO"), ("lumo_line", "LUMO")]:
        g.Add("xy", name=name)
        g.Set(f"{name}/xData", f"{name}_x")
        g.Set(f"{name}/yData", f"{name}_y")
        g.Set(f"{name}/marker", "none")
        g.Set(f"{name}/PlotLine/color", "darkgrey")
        g.Set(f"{name}/PlotLine/style", "dashed")
        g.Set(f"{name}/key", key)

    g.Add("label", name="gap_label",
          label=f"gap = {homo_lumo['gap_eV']:.2f} eV")
    g.Set("gap_label/xPos", 0.5)
    g.Set("gap_label/yPos", 0.92)
    g.Set("gap_label/alignHorz", "centre")


def add_multiline_graph(g, name, title, curves, xlabel, ylabel, emin, emax,
                         xdata="energy_eV", homo_lumo=None):
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

    if homo_lumo is not None:
        add_homo_lumo_lines(g, homo_lumo)

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
    parser.add_argument("--element4-label", default="O",
                         help="Label of the 4th projected element in dos_dataset.dat "
                              "(must match the pDOS_<label> column from parse_dos.py)")
    parser.add_argument("--homo-lumo", default=None,
                         help="Path to homo_lumo_dataset.dat (parse_homo_lumo.py). "
                              "If given, overlays HOMO/LUMO lines + gap label on the TDOS panel.")
    args = parser.parse_args()

    homo_lumo = read_homo_lumo(args.homo_lumo) if args.homo_lumo else None

    g = veusz.Embedded("dos_figure", hidden=True)

    g.To(g.Add("page", name="page1"))
    g.To(g.Add("grid", name="grid1", columns=2))

    g.ImportFile(args.data, "", linked=True)
    g.ImportFile(args.jdos, "", linked=True)

    add_multiline_graph(
        g, "graph_tdos", "TDOS",
        [("TDOS", "Total", "black"), ("DOS_s", "s", "red"), ("DOS_p", "p", "blue")],
        "Energy (eV)", "DOS (states/eV)", args.emin, args.emax,
        homo_lumo=homo_lumo
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
        [("pDOS_C", "C", "black"), ("pDOS_H", "H", "grey"),
         (f"pDOS_{args.element4_label}", args.element4_label, "red")],
        "Energy (eV)", "PDOS (states/eV)", args.emin, args.emax
    )

    add_placeholder_graph(g, "graph_reserved", "(reserved)")

    g.Save(f"{args.output}.vsz")
    g.Export(f"{args.output}.pdf")
    print(f"Saved -> {args.output}.vsz and {args.output}.pdf")
    g.Close()


if __name__ == "__main__":
    main()