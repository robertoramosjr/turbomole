#!/usr/bin/env python3
"""
count_wigner_window.py

Reports how many Wigner-sampled initial conditions fall inside a laser
energy window centered on a given excited state's ensemble band maximum
-- the practical question behind "how many initial conditions can I
actually excite into this band" (see parse_wigner_ic_dump.py for how the
input dataset is built from a split Newton-X initial-condition batch).

The band maximum is located by building an oscillator-strength-weighted,
Gaussian-broadened ensemble spectrum for the target state (same
convention as parse_excitations.py's JDOS) and taking its argmax --
NOT the single-point (Franck-Condon) vertical energy, since states like
azobenzene's n->pi* are symmetry-forbidden at equilibrium and all of the
ensemble's absorption intensity comes from distorted geometries.

For each window half-width and oscillator-strength cutoff, reports:
  - how many points fall inside the energy window at all
  - how many of those also clear the oscillator-strength cutoff
  - the distribution (min/25%/median/75%/max/mean) of oscillator
    strength among the in-window points, not just a pass/fail count

Usage:
    python ~/scripts/count_wigner_window.py \
        --dump wigner_1000_states_dump.csv --state 1 \
        --window-widths 0.05,0.10,0.20,0.30 \
        --osc-cutoffs 0,1e-4,1e-3 \
        --osc-column osc_len --broadening 0.10 --npoints 4000 \
        --output npi_star_window_counts.csv
"""

import argparse
import csv

import numpy as np


def read_dump(path, state, osc_column):
    energies, osc = [], []
    with open(path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if int(row["state"]) != state:
                continue
            energies.append(float(row["energy_eV"]))
            osc.append(float(row[osc_column]))
    if not energies:
        raise ValueError(f"No rows found for state={state} in {path}")
    return np.array(energies), np.array(osc)


def gaussian_broadened_spectrum(energies, osc, broadening_eV, npoints):
    emin = max(0.0, energies.min() - 5 * broadening_eV)
    emax = energies.max() + 5 * broadening_eV
    grid = np.linspace(emin, emax, npoints)
    profile = np.zeros_like(grid)
    for e0, f in zip(energies, osc):
        profile += f * np.exp(-((grid - e0) ** 2) / (2 * broadening_eV ** 2))
    return grid, profile


def summarize_osc(values):
    if len(values) == 0:
        return {"min": float("nan"), "p25": float("nan"), "median": float("nan"),
                "p75": float("nan"), "max": float("nan"), "mean": float("nan")}
    return {
        "min": float(np.min(values)), "p25": float(np.percentile(values, 25)),
        "median": float(np.median(values)), "p75": float(np.percentile(values, 75)),
        "max": float(np.max(values)), "mean": float(np.mean(values)),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Count Wigner initial conditions inside a laser energy window "
                    "centered on an excited state's ensemble band maximum.")
    parser.add_argument("--dump", required=True, help="Path to the parse_wigner_ic_dump.py CSV")
    parser.add_argument("--state", type=int, required=True,
                         help="Turbomole $soes state index to analyze (e.g. 1 for S1/n->pi*)")
    parser.add_argument("--window-widths", default="0.05,0.10,0.20,0.30",
                         help="Comma-separated half-widths in eV (window = center +/- width)")
    parser.add_argument("--osc-cutoffs", default="0,1e-4,1e-3",
                         help="Comma-separated oscillator-strength cutoffs (f > cutoff)")
    parser.add_argument("--osc-column", default="osc_len", choices=["osc_len", "osc_vel"],
                         help="Which oscillator-strength gauge to use")
    parser.add_argument("--broadening", type=float, default=0.10,
                         help="Gaussian broadening (eV) used only to locate the ensemble "
                              "band maximum, not for the window counting itself")
    parser.add_argument("--npoints", type=int, default=4000,
                         help="Grid resolution for the broadened spectrum used to find the peak")
    parser.add_argument("--output", default="wigner_window_counts.csv")
    args = parser.parse_args()

    window_widths = [float(w) for w in args.window_widths.split(",")]
    osc_cutoffs = [float(c) for c in args.osc_cutoffs.split(",")]

    energies, osc = read_dump(args.dump, args.state, args.osc_column)
    print(f"Loaded {len(energies)} points for state {args.state} "
          f"(energy range {energies.min():.3f}-{energies.max():.3f} eV, "
          f"{args.osc_column} range {osc.min():.2e}-{osc.max():.2e})")

    grid, profile = gaussian_broadened_spectrum(energies, osc, args.broadening, args.npoints)
    band_max_eV = grid[np.argmax(profile)]
    print(f"Ensemble band maximum (broadening={args.broadening} eV): {band_max_eV:.4f} eV")

    rows = []
    for width in window_widths:
        in_window = np.abs(energies - band_max_eV) <= width
        n_in_window = int(in_window.sum())
        osc_in_window = osc[in_window]
        stats = summarize_osc(osc_in_window)

        print(f"\nWindow +/-{width:.2f} eV around {band_max_eV:.4f} eV: "
              f"{n_in_window} points")
        print(f"  osc.strength in window: min={stats['min']:.2e} p25={stats['p25']:.2e} "
              f"median={stats['median']:.2e} p75={stats['p75']:.2e} max={stats['max']:.2e} "
              f"mean={stats['mean']:.2e}")

        for cutoff in osc_cutoffs:
            n_above_cutoff = int((osc_in_window > cutoff).sum())
            print(f"  f > {cutoff:.0e}: {n_above_cutoff} points")
            rows.append({
                "window_width_eV": width, "band_max_eV": band_max_eV,
                "n_in_window": n_in_window, "osc_cutoff": cutoff,
                "n_above_cutoff": n_above_cutoff,
                "osc_min": stats["min"], "osc_p25": stats["p25"], "osc_median": stats["median"],
                "osc_p75": stats["p75"], "osc_max": stats["max"], "osc_mean": stats["mean"],
            })

    fieldnames = ["window_width_eV", "band_max_eV", "n_in_window", "osc_cutoff",
                  "n_above_cutoff", "osc_min", "osc_p25", "osc_median", "osc_p75",
                  "osc_max", "osc_mean"]
    with open(args.output, "w") as f:
        f.write(",".join(fieldnames) + "\n")
        for r in rows:
            f.write(",".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k])
                              for k in fieldnames) + "\n")
    print(f"\nTable saved -> {args.output}")


if __name__ == "__main__":
    main()
