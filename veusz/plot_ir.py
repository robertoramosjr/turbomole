#!/usr/bin/env python3
"""
plot_ir.py

Plots the simulated IR spectrum from Turbomole aoforce output:
Lorentzian-broadened envelope (continuous curve, matching the visual
style of experimental/simulated IR spectra like Fig. 7 of Januario &
Cabral, 2026) overlaid with the raw stick spectrum (discrete vibrational
modes). Broadening width is an assumption, not from the paper (which
doesn't report its exact IR broadening parameter) -- flagged here.

Usage:
    python ~/scripts/veusz/plot_ir.py --data ir_dataset.dat \
        --output ir_figure --fwhm 10 --npoints 3000
"""

import argparse
import os

import numpy as np

if not os.environ.get("DISPLAY"):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import veusz.embed as veusz


def lorentzian_broaden(wavenumbers, intensities, fwhm, npoints, wmin=None, wmax=None):
    gamma = fwhm / 2.0
    if wmin is None:
        wmin = max(0.0, wavenumbers.min() - 10 * fwhm)
    if wmax is None:
        wmax = wavenumbers.max() + 10 * fwhm

    grid = np.linspace(wmin, wmax, npoints)
    profile = np.zeros_like(grid)
    for w0, inten in zip(wavenumbers, intensities):
        profile += inten * (gamma ** 2) / ((grid - w0) ** 2 + gamma ** 2)

    return grid, profile


def main():
    parser = argparse.ArgumentParser(description="Plot the TD-DFT/aoforce IR spectrum in Veusz.")
    parser.add_argument("--data", required=True, help="Path to ir_dataset.dat")
    parser.add_argument("--output", default="ir_figure", help="Output filename prefix")
    parser.add_argument("--fwhm", type=float, default=10.0,
                         help="Lorentzian FWHM (cm-1) for the broadened envelope")
    parser.add_argument("--npoints", type=int, default=3000)
    args = parser.parse_args()

    raw = np.loadtxt(args.data, skiprows=1)
    modes, wavenumbers, intensities = raw[:, 0], raw[:, 1], raw[:, 2]

    grid, profile = lorentzian_broaden(wavenumbers, intensities, args.fwhm, args.npoints)

    broadened_path = f"{args.output}_broadened.dat"
    with open(broadened_path, "w") as f:
        f.write("descriptor ir_wavenumber_cm1,ir_broadened_intensity\n")
        for w, i in zip(grid, profile):
            f.write(f"{w:.4f} {i:.6f}\n")

    g = veusz.Embedded("ir_figure", hidden=True)

    g.To(g.Add("page", name="page1"))
    g.To(g.Add("graph", name="graph1"))

    g.ImportFile(args.data, "", linked=True)
    g.ImportFile(broadened_path, "", linked=True)

    g.Set("x/label", "Wavenumber (cm-1)")
    g.Set("x/min", float(grid.min()))
    g.Set("x/max", float(grid.max()))

    g.Set("y/label", "IR intensity (km/mol)")
    g.Set("y/min", 0)

    # Broadened envelope
    g.Add("xy", name="xy_envelope")
    g.Set("xy_envelope/xData", "ir_wavenumber_cm1")
    g.Set("xy_envelope/yData", "ir_broadened_intensity")
    g.Set("xy_envelope/marker", "none")
    g.Set("xy_envelope/PlotLine/color", "black")
    g.Set("xy_envelope/key", "Broadened envelope")

    # Raw stick spectrum as thin bars
    g.Add("bar", name="bar_sticks")
    g.Set("bar_sticks/lengths", ["ir_intensity_kmmol"])
    g.Set("bar_sticks/posn", "wavenumber_cm1")
    g.Set("bar_sticks/BarFill/fills", [("solid", "red", False)])
    g.Set("bar_sticks/keys", ["Stick spectrum"])

    g.Add("key", name="key1")

    g.Save(f"{args.output}.vsz")
    g.Export(f"{args.output}.pdf")
    print(f"Saved -> {args.output}.vsz and {args.output}.pdf")
    g.Close()


if __name__ == "__main__":
    main()