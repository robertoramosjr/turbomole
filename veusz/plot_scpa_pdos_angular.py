#!/usr/bin/env python3
"""
plot_scpa_pdos_angular.py

Plots the SCPA-partitioned, angular-momentum-resolved PDOS (s/p per
element: C, O, H) produced by parse_scpa_pdos_angular.py, all overlaid
in a single graph -- color by element (same palette as
plot_becke_pdos.py: O red, C blue, H green), line style by angular
momentum (solid = s, dashed = p). TDOS is shown in black for reference.

Note: this fragment definition only resolves s/p functions -- any d (or
higher) polarization functions on C/O aren't assigned to any of the 6
tracked fragments, so pDOS_Cs+pDOS_Cp+pDOS_Os+pDOS_Op+pDOS_Hs+pDOS_Hp
does not sum exactly to TDOS; the gap is the d-shell contribution.

Whatever Gaussian FWHM was used when this DOS_curve.txt was exported
from Multiwfn controls how smooth the curves already are -- if it was
narrow (individual orbitals resolved as near-delta spikes on the energy
grid, rather than a smooth envelope), pass --sigma to additionally
broaden on top of that at plot time, same technique as
plot_qpenergies_dos.py.

CAUTION: keep --sigma well below the HOMO-LUMO gap width, or the smeared
tail of the frontier peaks bleeds visibly into the gap and makes the
HOMO/LUMO lines look like they're cutting through occupied density
instead of sitting exactly on the terminal peaks (confirmed on the
artepillin_C export: --sigma 0.4 did this against a 3.68 eV gap;
--sigma 0.1 kept the frontier peaks sharp and isolated).

Usage:
    python ~/scripts/veusz/plot_scpa_pdos_angular.py \
        --data dos_scpa_angular_dataset.dat --output scpa_pdos_figure

    # extra Gaussian smearing on top of whatever's already in the data
    # (0 disables; 0.1 eV made a very spiky artepillin_C export readable
    # without smearing the frontier peaks into the gap -- start there):
    python ~/scripts/veusz/plot_scpa_pdos_angular.py \
        --data dos_scpa_angular_dataset.dat --sigma 0.1 --output scpa_pdos_figure

    # optionally mark HOMO/LUMO and annotate the gap, same convention
    # as plot_becke_pdos.py (derived from Multiwfn's orginfo.txt):
    python ~/scripts/veusz/plot_scpa_pdos_angular.py \
        --data dos_scpa_angular_dataset.dat --orginfo orginfo.txt \
        --output scpa_pdos_figure
"""

import argparse
import os

import numpy as np

if not os.environ.get("DISPLAY"):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import veusz.embed as veusz

HARTREE_TO_EV = 27.211386245988


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
        description="Plot SCPA angular-momentum-resolved PDOS (s/p per "
                    "element) in a single overlaid graph.")
    parser.add_argument("--data", required=True, help="Path to dos_scpa_angular_dataset.dat")
    parser.add_argument("--output", default="scpa_pdos_figure", help="Output filename prefix")
    parser.add_argument("--emin", type=float, default=-20.0)
    parser.add_argument("--emax", type=float, default=10.0)
    parser.add_argument("--orginfo",
                         help="Path to Multiwfn's orginfo.txt, used to derive HOMO/LUMO "
                              "automatically (optional, unless --homo/--lumo are given)")
    parser.add_argument("--homo", type=float,
                         help="HOMO energy in eV -- overrides --orginfo if given")
    parser.add_argument("--lumo", type=float,
                         help="LUMO energy in eV -- overrides --orginfo if given")
    parser.add_argument("--sigma", type=float, default=0.0,
                         help="Standard deviation (eV) of an additional Gaussian "
                              "smearing applied to TDOS and all six PDOS curves "
                              "before plotting, on top of whatever broadening is "
                              "already baked into --data (0, the default, leaves "
                              "the data untouched)")
    args = parser.parse_args()

    homo = lumo = None
    if args.homo is not None and args.lumo is not None:
        homo, lumo = args.homo, args.lumo
    elif args.orginfo:
        homo, lumo = homo_lumo_from_orginfo(args.orginfo)
        print(f"HOMO/LUMO derived from {args.orginfo}")

    labels = ["Cs", "Cp", "Os", "Op", "Hs", "Hp"]
    columns = np.loadtxt(args.data, skiprows=1, unpack=True)
    energy, tdos = columns[0], columns[1]
    pdos = dict(zip(labels, columns[2:8]))

    data_path = args.data
    if args.sigma > 0:
        dx = float(energy[1] - energy[0])
        tdos = gaussian_smear(tdos, dx, args.sigma)
        pdos = {name: gaussian_smear(values, dx, args.sigma)
                for name, values in pdos.items()}

        data_path = f"{args.output}_smeared.dat"
        with open(data_path, "w") as f:
            f.write("descriptor energy_eV,TDOS," +
                    ",".join(f"pDOS_{name}" for name in labels) + "\n")
            for i in range(len(energy)):
                row = [f"{energy[i]:.6f}", f"{tdos[i]:.6f}"] + \
                      [f"{pdos[name][i]:.6f}" for name in labels]
                f.write(" ".join(row) + "\n")

    in_window = (energy >= args.emin) & (energy <= args.emax)
    ymax = float(max(tdos[in_window].max(),
                      max(v[in_window].max() for v in pdos.values())))

    g = veusz.Embedded("scpa_pdos_figure", hidden=True)
    g.EnableToolbar(False)

    page = g.Root.Add("page", name="page1")
    graph = page.Add("graph", name="graph1")

    g.To(graph.path)
    g.ImportFile(data_path, "", linked=True)

    # Set properties on the graph's default x/y axes rather than adding new
    # ones -- graph.Add("axis", ...) creates a duplicate ghost axis on top
    # of the auto-created default, doubling up the tick labels.
    graph.x.label.val = "Energy (eV)"
    graph.x.min.val = args.emin
    graph.x.max.val = args.emax

    graph.y.label.val = "DOS (arb. units)"
    graph.y.min.val = 0
    if homo is not None:
        # Explicit headroom above the tallest curve so the HOMO/LUMO/gap
        # labels (placed in axis coordinates) aren't clipped by an
        # auto-ranged top.
        graph.y.max.val = ymax * 1.2

    curves = [
        ("xy_tdos", "TDOS", "Total", "black", "solid", "2pt"),
        ("xy_cs", "pDOS_Cs", "C (s)", "blue", "solid", "1pt"),
        ("xy_cp", "pDOS_Cp", "C (p)", "blue", "dashed", "1pt"),
        ("xy_os", "pDOS_Os", "O (s)", "red", "solid", "1pt"),
        ("xy_op", "pDOS_Op", "O (p)", "red", "dashed", "1pt"),
        ("xy_hs", "pDOS_Hs", "H (s)", "green", "solid", "1pt"),
        ("xy_hp", "pDOS_Hp", "H (p)", "green", "dashed", "1pt"),
    ]
    for name, ydata, key, color, style, width in curves:
        xy = graph.Add("xy", name=name)
        xy.xData.val = "energy_eV"
        xy.yData.val = ydata
        xy.marker.val = "none"
        xy.PlotLine.color.val = color
        xy.PlotLine.style.val = style
        xy.PlotLine.width.val = width
        xy.key.val = key

    if homo is not None:
        gap = lumo - homo

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

        print(f"HOMO = {homo:.4f} eV, LUMO = {lumo:.4f} eV, gap = {gap:.4f} eV")

    graph.Add("key", name="key1")
    g.Set("key1/Text/size", "8pt")

    g.Save(f"{args.output}.vsz")
    g.Export(f"{args.output}.pdf")
    print(f"Saved -> {args.output}.vsz and {args.output}.pdf")
    g.Close()


if __name__ == "__main__":
    main()
