#!/usr/bin/env python3
"""
plot_qpenergies_dos.py

Stacks the Kohn-Sham (HSE06) and G0W0-corrected orbital density
datasets produced by parse_qpenergies_dos.py as two vertically aligned
panels sharing the same energy axis, for direct visual comparison of
the frontier-region gap opening under the GW quasiparticle correction.

Usage:
    python ~/scripts/veusz/plot_qpenergies_dos.py --ks dos_ks_dataset.dat \
        --qp dos_qp_dataset.dat --output qpdos_figure --emin -20 --emax 10 \
        --sigma 0.01
"""

import argparse
import os

import numpy as np

if not os.environ.get("DISPLAY"):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import veusz.embed as veusz


def gaussian_smear(y, dx, sigma):
    """Convolve y (sampled on a uniform grid with spacing dx) with a
    Gaussian kernel of the given standard deviation (same units as dx)."""
    if sigma <= 0:
        return y
    half_width = max(1, int(np.ceil(5 * sigma / dx)))
    xs = np.arange(-half_width, half_width + 1) * dx
    kernel = np.exp(-(xs ** 2) / (2 * sigma ** 2))
    kernel /= kernel.sum()
    return np.convolve(y, kernel, mode="same")


def write_smeared_dataset(path, xname, xvalues, yname, yvalues):
    with open(path, "w") as f:
        f.write(f"descriptor {xname},{yname}\n")
        for x, y in zip(xvalues, yvalues):
            f.write(f"{x:.6f} {y:.6f}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Plot the KS (HSE06) vs. G0W0-corrected orbital density overlay in Veusz.")
    parser.add_argument("--ks", required=True, help="Path to dos_ks_dataset.dat")
    parser.add_argument("--qp", required=True, help="Path to dos_qp_dataset.dat")
    parser.add_argument("--output", default="qpdos_figure", help="Output filename prefix")
    parser.add_argument("--emin", type=float, default=-20.0)
    parser.add_argument("--emax", type=float, default=10.0)
    parser.add_argument("--sigma", type=float, default=0.01,
                         help="Standard deviation (eV) of an additional Gaussian "
                              "smearing applied to the KS/QP density curves before "
                              "plotting, on top of whatever broadening is already "
                              "baked into the input datasets (0 disables)")
    args = parser.parse_args()

    ks_energy, ks_density = np.loadtxt(args.ks, skiprows=1, unpack=True)
    qp_energy, qp_density = np.loadtxt(args.qp, skiprows=1, unpack=True)

    ks_density = gaussian_smear(ks_density, ks_energy[1] - ks_energy[0], args.sigma)
    qp_density = gaussian_smear(qp_density, qp_energy[1] - qp_energy[0], args.sigma)

    ks_smeared_path = f"{args.output}_ks_smeared.dat"
    write_smeared_dataset(ks_smeared_path, "energy_eV", ks_energy, "ks_density", ks_density)

    qp_smeared_path = f"{args.output}_qp_smeared.dat"
    write_smeared_dataset(qp_smeared_path, "energy_eV_qp", qp_energy, "qp_density", qp_density)

    g = veusz.Embedded("qpdos_figure", hidden=True)
    g.EnableToolbar(False)

    page = g.Root.Add("page", name="page1")
    grid = page.Add("grid", name="grid1", columns=1)

    g.To(grid.path)
    g.ImportFile(ks_smeared_path, "", linked=True)
    g.ImportFile(qp_smeared_path, "", linked=True)

    # Top panel: Kohn-Sham (HSE06)
    graph_ks = grid.Add("graph", name="graph_ks")
    g.To(graph_ks.path)

    # Set properties on the graph's default x/y axes rather than adding new
    # ones -- graph.Add("axis", ...) creates a duplicate ghost axis on top
    # of the auto-created default, doubling up the tick labels.
    graph_ks.x.label.val = "Energy (eV)"
    graph_ks.x.min.val = args.emin
    graph_ks.x.max.val = args.emax

    graph_ks.y.label.val = "KS density (arb. units)"
    graph_ks.y.min.val = "Auto"

    xy_ks = graph_ks.Add("xy", name="xy_ks")
    xy_ks.xData.val = "energy_eV"
    xy_ks.yData.val = "ks_density"
    xy_ks.marker.val = "none"
    xy_ks.PlotLine.color.val = "black"
    xy_ks.PlotLine.width.val = "1.5pt"

    label_ks = graph_ks.Add("label", name="title_label")
    label_ks.label.val = "Kohn-Sham (HSE06)"
    label_ks.xPos.val = 0.97
    label_ks.yPos.val = 0.9
    label_ks.alignHorz.val = "right"

    # Bottom panel: G0W0-corrected (QP)
    graph_qp = grid.Add("graph", name="graph_qp")
    g.To(graph_qp.path)

    graph_qp.x.label.val = "Energy (eV)"
    graph_qp.x.min.val = args.emin
    graph_qp.x.max.val = args.emax

    graph_qp.y.label.val = "QP density (arb. units)"
    graph_qp.y.min.val = "Auto"

    xy_qp = graph_qp.Add("xy", name="xy_qp")
    xy_qp.xData.val = "energy_eV_qp"
    xy_qp.yData.val = "qp_density"
    xy_qp.marker.val = "none"
    xy_qp.PlotLine.color.val = "red"
    xy_qp.PlotLine.width.val = "1.5pt"

    label_qp = graph_qp.Add("label", name="title_label")
    label_qp.label.val = "G0W0-corrected (QP)"
    label_qp.xPos.val = 0.97
    label_qp.yPos.val = 0.9
    label_qp.alignHorz.val = "right"

    g.Save(f"{args.output}.vsz")
    g.Export(f"{args.output}.pdf")
    print(f"Saved -> {args.output}.vsz and {args.output}.pdf")
    g.Close()


if __name__ == "__main__":
    main()