#!/usr/bin/env python3
"""
build_population_weighted_spectrum.py

Layer A's main deliverable: combines the per-conformer IR
(vibspectrum, from aoforce) or UV-Vis (escf.out, TD-DFT) spectra of
every population-filtered conformer of one protonation state into a
single population-weighted broadened envelope, using the r2SCAN-3c
population weights from rerank_layer_a_dft.py's output
(population_table_dft.csv).

Each conformer's stick spectrum is broadened individually (same
formulas/defaults as plot_ir.py / parse_excitations.py: Lorentzian for
IR, Gaussian for UV-Vis), scaled by that conformer's population
weight, and summed on a common grid -- this is the discretized version
of "average the spectrum over the Boltzmann ensemble."

IR-specific caveat (read before trusting the output): the geometries
were optimized at r2SCAN-3c/def2-mTZVP with a loosened, hard-capped
convergence (not a true stationary point even on that surface), and
the Hessian is computed at a different level (PBE0/def2-SVP) --
confirmed empirically that this leaves EVERY conformer with several
genuine imaginary modes (large negative wavenumbers, not just
near-zero numerical noise), not just a few difficult cases. Real
translation/rotation modes come out essentially exactly at 0 cm-1 in
this data; imaginary modes are the substantially negative ones (seen:
down to -273 cm-1) that sort *before* the near-zero translation/
rotation modes in Turbomole's ascending-frequency ordering -- so
`parse_ir.py`'s "exclude the first 6 by position" convention does NOT
correctly separate them here. This script instead filters by value: a
mode is kept only if its wavenumber exceeds --min-wavenumber (default
10 cm-1, chosen from this dataset: translation/rotation modes sit at
~0.00 cm-1, imaginary modes are substantially negative, and the lowest
real vibrational mode observed was ~41 cm-1 -- adjust for other
molecules if their lowest genuine vibrational mode is lower than that).
Imaginary modes are dropped, not modeled as negative-frequency bands
(they don't correspond to any real absorption) -- this is a documented
approximation of the screening stage, not the final production
spectrum (see protonation_and_conformer_investigation.md, Estagio 8).

Usage:
    python ~/work_turbomole/scripts/build_population_weighted_spectrum.py \
        --kind ir \
        --population-table layer_a/population_table_dft.csv \
        --jobs-dir layer_a/pbe0_jobs \
        --output layer_a/ir_population_weighted.dat

    python ~/work_turbomole/scripts/build_population_weighted_spectrum.py \
        --kind uvvis \
        --population-table layer_a/population_table_dft.csv \
        --jobs-dir layer_a/pbe0_jobs \
        --output layer_a/uvvis_population_weighted.dat
"""

import argparse
import csv
import os
import re

import numpy as np

VIB_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(?:([a-zA-Z0-9']+)\s+)?([-\d.]+)\s+([-\d.]+)\s+(\S+)\s+(\S+)\s*$"
)
ESCF_ROW_RE = re.compile(
    r"^\s*\|\s*(\d+)\s+a\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
    r"\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
)


def read_population_table(path):
    """Return list of (job_dir, population) for kept conformers."""
    entries = []
    with open(path) as f:
        for row in csv.DictReader(f):
            entries.append((row["job_dir"], float(row["r2scan3c_population"])))
    return entries


def parse_ir_modes(vibspectrum_path, min_wavenumber):
    modes = []
    with open(vibspectrum_path) as f:
        for line in f:
            if line.strip().startswith(("$", "#")):
                continue
            m = VIB_ROW_RE.match(line)
            if not m:
                continue
            wavenumber = float(m.group(3))
            intensity = float(m.group(4))
            if wavenumber > min_wavenumber:
                modes.append((wavenumber, intensity))
    return modes


def parse_uvvis_states(escf_path):
    states = []
    with open(escf_path) as f:
        for line in f:
            m = ESCF_ROW_RE.match(line)
            if not m:
                continue
            energy_eV = float(m.group(2))
            osc_len = float(m.group(6))
            states.append((energy_eV, osc_len))
    return states


def lorentzian_broaden(sticks, fwhm, npoints, xmin, xmax):
    gamma = fwhm / 2.0
    grid = np.linspace(xmin, xmax, npoints)
    profile = np.zeros_like(grid)
    for x0, height in sticks:
        profile += height * (gamma ** 2) / ((grid - x0) ** 2 + gamma ** 2)
    return grid, profile


def gaussian_broaden(sticks, sigma, npoints, xmin, xmax):
    grid = np.linspace(xmin, xmax, npoints)
    profile = np.zeros_like(grid)
    for x0, height in sticks:
        profile += height * np.exp(-((grid - x0) ** 2) / (2 * sigma ** 2))
    return grid, profile


def main():
    parser = argparse.ArgumentParser(
        description="Combine per-conformer IR or UV-Vis spectra into one "
                    "population-weighted broadened envelope for a protonation state.")
    parser.add_argument("--kind", required=True, choices=["ir", "uvvis"])
    parser.add_argument("--population-table", required=True,
                         help="rerank_layer_a_dft.py output "
                              "(population_table_dft.csv)")
    parser.add_argument("--jobs-dir", required=True,
                         help="Parent folder of job_NNNN/ dirs with "
                              "vibspectrum or escf.out (pbe0_jobs)")
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-wavenumber", type=float, default=10.0,
                         help="[ir only] minimum wavenumber (cm-1) to keep as "
                              "a real vibrational mode -- excludes both "
                              "near-zero translation/rotation and imaginary "
                              "modes [default: 10.0]")
    parser.add_argument("--fwhm", type=float, default=10.0,
                         help="[ir only] Lorentzian FWHM in cm-1 [default: 10.0, "
                              "matches plot_ir.py]")
    parser.add_argument("--broadening", type=float, default=0.4,
                         help="[uvvis only] Gaussian sigma in eV [default: 0.4, "
                              "matches parse_excitations.py]")
    parser.add_argument("--npoints", type=int, default=3000)
    parser.add_argument("--xmin", type=float, default=None)
    parser.add_argument("--xmax", type=float, default=None)
    args = parser.parse_args()

    entries = read_population_table(args.population_table)
    total_weight = sum(w for _, w in entries)
    print(f"Loaded {len(entries)} conformers, total population = {total_weight:.4f} "
          f"(should be close to the Layer A cumulative cutoff, e.g. ~0.95)")

    all_sticks_by_conf = []
    dropped_imaginary_count = 0
    for job_dir, weight in entries:
        job_path = os.path.join(args.jobs_dir, job_dir)
        if args.kind == "ir":
            vib_path = os.path.join(job_path, "vibspectrum")
            modes = parse_ir_modes(vib_path, args.min_wavenumber)
            with open(vib_path) as f:
                n_imaginary = sum(
                    1 for line in f
                    if (m := VIB_ROW_RE.match(line)) and float(m.group(3)) < 0
                )
            dropped_imaginary_count += n_imaginary
            all_sticks_by_conf.append((weight, modes))
        else:
            escf_path = os.path.join(job_path, "escf.out")
            states = parse_uvvis_states(escf_path)
            all_sticks_by_conf.append((weight, states))

    if args.kind == "ir":
        print(f"Dropped {dropped_imaginary_count} imaginary modes total across "
              f"the ensemble (not modeled as bands -- see script docstring)")

    all_x = [x for _, sticks in all_sticks_by_conf for x, _ in sticks]
    if not all_x:
        raise SystemExit("No modes/states survived parsing -- nothing to plot.")
    pad = args.fwhm * 10 if args.kind == "ir" else args.broadening * 5
    xmin = args.xmin if args.xmin is not None else max(0.0, min(all_x) - pad)
    xmax = args.xmax if args.xmax is not None else max(all_x) + pad

    combined = None
    grid = None
    for weight, sticks in all_sticks_by_conf:
        if not sticks:
            continue
        if args.kind == "ir":
            grid, profile = lorentzian_broaden(sticks, args.fwhm, args.npoints, xmin, xmax)
        else:
            grid, profile = gaussian_broaden(sticks, args.broadening, args.npoints, xmin, xmax)
        weighted = weight * profile
        combined = weighted if combined is None else combined + weighted

    xlabel = "wavenumber_cm1" if args.kind == "ir" else "energy_eV"
    ylabel = "ir_intensity_weighted" if args.kind == "ir" else "absorption_weighted"
    with open(args.output, "w") as f:
        f.write(f"descriptor {xlabel},{ylabel}\n")
        for x, y in zip(grid, combined):
            f.write(f"{x:.6f} {y:.8f}\n")

    print(f"Population-weighted {args.kind.upper()} envelope written -> {args.output}")


if __name__ == "__main__":
    main()
