#!/usr/bin/env python3
"""
predict_uvvis_vs_ph.py

Predicted UV-Vis spectrum of artepilin C as a function of pH
(protonation_and_conformer_investigation.md Sec. 9.37): the
population-weighted HSE06/def2-TZVP TD-DFT spectrum of each protonation
microstate (REAL 60-state escf.out sticks of every conformer that was
computed -- no surrogate, no S1-only truncation), mixed with the
speciation fractions from plot_speciation_vs_ph.py:

    S(E; pH) = sum_i f_i(pH) * S_i(E)

Per-state spectra: Gaussian-broadened sticks (sigma, default 0.4 eV --
the project convention of build_population_weighted_spectrum.py; note
this is broad compared with a typical experimental band), weights =
r2SCAN-3c populations renormalized WITHIN the set of conformers that
have a real escf.out (neutral 37/170 = 50% of the population,
monoanion A 20/143 = 27%, monoanion B 11/11, dianion 3/3). The
coverage uncertainty of that renormalization is quantified in
estimate_model_gap_impact.py (a few 0.01 eV for the robust observables
of monoanion A).

Reports, per state and per pH, the band MAXIMUM and the f-weighted
CENTROID over --window (default 3.0-5.5 eV). This is the prediction to
be compared with the experimental near-UV band position vs pH
(Camuri/Pazin titration); the experimental spectra are NOT available on
this machine, so no comparison is made here.

Usage:
    python ~/work_turbomole/scripts/predict_uvvis_vs_ph.py \
        --project-root . --speciation speciation_v3_4refs.dat \
        --ph 2 4.65 7.4 10 12 --output-prefix uvvis_vs_ph
"""

import argparse
import csv
import os
import re

import numpy as np
import matplotlib.pyplot as plt

ESCF_ROW_RE = re.compile(
    r"^\s*\|\s*(\d+)\s+a\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
    r"\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
)
STATE_DIRS = {"N": "neutral", "A": "monoanion_a", "B": "monoanion_b", "D": "dianion"}
COLORS = {"N": "#1f77b4", "A": "#2ca02c", "B": "#d62728", "D": "#9467bd"}


def sticks(path):
    st = []
    with open(path) as f:
        for line in f:
            m = ESCF_ROW_RE.match(line)
            if m:
                st.append((float(m.group(3)), float(m.group(7))))  # eV, f(length)
    return st


def state_spectrum(root, state, grid, sigma):
    d = os.path.join(root, f"artepillin_c_{state}_d3_hse06")
    with open(os.path.join(d, "layer_a", "population_table_dft.csv")) as f:
        pop = {r["job_dir"]: float(r["r2scan3c_population"]) for r in csv.DictReader(f)}
    found = {}
    for sub in ("labeled_jobs", "confirmation_batch_jobs", "random_control_batch_jobs"):
        base = os.path.join(d, "layer_b", sub)
        if os.path.isdir(base):
            for j in os.listdir(base):
                p = os.path.join(base, j, "escf.out")
                if os.path.isfile(p) and j not in found and j in pop:
                    found[j] = p
    tot = sum(pop[j] for j in found)
    prof = np.zeros_like(grid)
    for j, p in found.items():
        w = pop[j] / tot
        for e, f in sticks(p):
            prof += w * f * np.exp(-((grid - e) ** 2) / (2 * sigma ** 2))
    cov = tot / sum(pop.values())
    return prof, len(found), len(pop), cov


def main():
    parser = argparse.ArgumentParser(
        description="Predicted population-weighted UV-Vis spectrum vs pH "
                    "(HSE06 real sticks x speciation fractions).")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--speciation", required=True, help="plot_speciation_vs_ph.py .dat")
    parser.add_argument("--ph", type=float, nargs="+", default=[2, 4.65, 7.4, 10, 12])
    parser.add_argument("--sigma", type=float, default=0.4)
    parser.add_argument("--window", type=float, nargs=2, default=[3.0, 5.5])
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    grid = np.linspace(1.5, 7.0, 2200)
    spec = {}
    print(f"{'state':5s} {'N conf':>10s} {'pop cov':>8s} {'max (eV/nm)':>16s} {'centroid (eV/nm)':>20s}")
    win = (grid >= args.window[0]) & (grid <= args.window[1])
    for k, s in STATE_DIRS.items():
        prof, n, ntot, cov = state_spectrum(args.project_root, s, grid, args.sigma)
        spec[k] = prof
        wmax = grid[win][np.argmax(prof[win])]
        cen = np.sum(grid[win] * prof[win]) / np.sum(prof[win])
        print(f"{k:5s} {n:4d}/{ntot:<4d} {cov:8.3f} {wmax:8.3f}/{1239.84 / wmax:5.0f} {cen:12.3f}/{1239.84 / cen:5.0f}")

    frac = np.loadtxt(args.speciation, skiprows=1)
    ph_grid = frac[:, 0]
    fr = {k: frac[:, i + 1] for i, k in enumerate("NABD")}
    rows = []
    print(f"\n{'pH':>6s} {'f_N':>6s} {'f_A':>6s} {'f_D':>6s}  {'max (eV/nm)':>16s} {'centroid (eV/nm)':>20s}")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for k in "NABD":
        axes[0].plot(1239.84 / grid, spec[k], color=COLORS[k], label=k)
    axes[0].set_xlim(220, 500); axes[0].set_xlabel("Wavelength (nm)"); axes[0].set_ylabel("Population-weighted absorption (arb.)")
    axes[0].set_title("Per-microstate (HSE06/def2-TZVP, real sticks)"); axes[0].legend(title="state")
    for ph in args.ph:
        i = int(np.argmin(np.abs(ph_grid - ph)))
        mix = sum(fr[k][i] * spec[k] for k in "NABD")
        wmax = grid[win][np.argmax(mix[win])]
        cen = np.sum(grid[win] * mix[win]) / np.sum(mix[win])
        rows.append((ph, wmax, cen))
        print(f"{ph:6.2f} {fr['N'][i]:6.3f} {fr['A'][i]:6.3f} {fr['D'][i]:6.3f}  {wmax:8.3f}/{1239.84 / wmax:5.0f} {cen:12.3f}/{1239.84 / cen:5.0f}")
        axes[1].plot(1239.84 / grid, mix, label=f"pH {ph:g}")
    axes[1].set_xlim(220, 500); axes[1].set_xlabel("Wavelength (nm)"); axes[1].set_title("Speciation-weighted spectrum"); axes[1].legend()
    cen_curve = []
    for i in range(0, len(ph_grid), 5):
        mix = sum(fr[k][i] * spec[k] for k in "NABD")
        cen_curve.append((ph_grid[i], np.sum(grid[win] * mix[win]) / np.sum(mix[win])))
    cc = np.array(cen_curve)
    axes[2].plot(cc[:, 0], 1239.84 / cc[:, 1], color="k")
    axes[2].set_xlabel("pH"); axes[2].set_ylabel("f-weighted centroid, %g-%g eV (nm)" % tuple(args.window))
    axes[2].set_title("Predicted band position vs pH")
    fig.tight_layout(); fig.savefig(f"{args.output_prefix}.png", dpi=180)
    with open(f"{args.output_prefix}.dat", "w") as f:
        f.write("descriptor pH,band_max_eV,centroid_eV\n")
        for ph, m_, c_ in rows:
            f.write(f"{ph:.2f} {m_:.4f} {c_:.4f}\n")
    print(f"\nWrote {args.output_prefix}.png/.dat")


if __name__ == "__main__":
    main()
