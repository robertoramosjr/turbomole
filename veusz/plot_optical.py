#!/usr/bin/env python3
"""
plot_optical.py

Plots the TD-DFT simulated absorption spectrum: Gaussian-broadened
envelope (continuous curve) overlaid with the raw stick spectrum.
Sticks are rendered as a single xy line with NaN-separated vertical
segments (one per excited state) -- the standard robust technique for
stick spectra in Veusz, avoiding uncertain bar-widget properties.

Usage:
    python ~/scripts/veusz/plot_optical.py --sticks excitation_sticks.dat \
        --jdos jdos_dataset.dat --output optical_figure
"""

import argparse
import numpy as np
import veusz.embed as veusz


def build_stick_line_dataset(sticks_path, output_path):
    """
    Reads excitation_sticks.dat (state, energy_eV, osc_len, osc_vel)
    and writes a NaN-separated line dataset: for each state, three
    rows (energy,0) -> (energy,osc_len) -> (NaN,NaN), so a single xy
    line plot renders isolated vertical sticks.
    """
    raw = np.loadtxt(sticks_path, skiprows=1)
    energies = raw[:, 1]
    osc_lens = raw[:, 2]

    with open(output_path, "w") as f:
        f.write("descriptor stick_x,stick_y\n")
        for e, o in zip(energies, osc_lens):
            f.write(f"{e:.6f} 0.0\n")
            f.write(f"{e:.6f} {o:.6f}\n")
            f.write("nan nan\n")


def main():
    parser = argparse.ArgumentParser(description="Plot TD-DFT absorption spectrum in Veusz.")
    parser.add_argument("--sticks", required=True, help="Path to excitation_sticks.dat")
    parser.add_argument("--jdos", required=True, help="Path to jdos_dataset.dat (broadened envelope)")
    parser.add_argument("--output", default="optical_figure", help="Output filename prefix")
    args = parser.parse_args()

    stick_line_path = f"{args.output}_sticks_line.dat"
    build_stick_line_dataset(args.sticks, stick_line_path)

    g = veusz.Embedded("optical_figure", hidden=True)

    g.To(g.Add("page", name="page1"))
    g.Set("width", "16cm")
    g.Set("height", "10cm")
    g.To(g.Add("graph", name="graph1"))

    g.ImportFile(stick_line_path, "", linked=True)
    g.ImportFile(args.jdos, "", linked=True)

    g.Set("x/label", "Energy (eV)")
    g.Set("x/min", 0)
    g.Set("x/max", 10)

    g.Set("y/label", "Intensity (arb. units)")

    # Stick spectrum: NaN-separated vertical line segments
    g.Add("xy", name="xy_sticks")
    g.Set("xy_sticks/xData", "stick_x")
    g.Set("xy_sticks/yData", "stick_y")
    g.Set("xy_sticks/marker", "none")
    g.Set("xy_sticks/PlotLine/color", "red")
    g.Set("xy_sticks/PlotLine/width", "1pt")
    g.Set("xy_sticks/key", "Stick spectrum (osc. strength)")

    # Broadened envelope
    g.Add("xy", name="xy_envelope")
    g.Set("xy_envelope/xData", "jdos_energy_eV")
    g.Set("xy_envelope/yData", "jdos_intensity")
    g.Set("xy_envelope/marker", "none")
    g.Set("xy_envelope/PlotLine/color", "black")
    g.Set("xy_envelope/PlotLine/width", "1.5pt")
    g.Set("xy_envelope/key", "Broadened envelope")

    g.Add("key", name="key1")

    g.Save(f"{args.output}.vsz")
    g.Export(f"{args.output}.pdf")
    print(f"Saved -> {args.output}.vsz and {args.output}.pdf")
    g.Close()


if __name__ == "__main__":
    main()