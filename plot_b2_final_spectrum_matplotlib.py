#!/usr/bin/env python3
"""
plot_b2_final_spectrum_matplotlib.py

Immediate-viewing fallback for the final B2 population-weighted UV-Vis
figure (neutral + monoanion_b), since Veusz is not installed anywhere
on this machine (checked all 7 conda envs, 2026-09-11) -- the existing
.vsz files in this project were generated on a different machine. This
produces the same plot with matplotlib so there is a real image to
look at now; plot_b2_final_spectrum_veusz.py produces the polished
Veusz version to run wherever Veusz is available.

Usage:
    python ~/work_turbomole/scripts/plot_b2_final_spectrum_matplotlib.py \
        --neutral-envelope b2_final_spectrum_neutral.dat \
        --monoanion-b-envelope b2_final_spectrum_monoanion_b.dat \
        --mismatch-rate 0.083 --mismatch-ci-low 0.0287 --mismatch-ci-high 0.2183 \
        --output b2_final_spectrum.png
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt


def load_envelope(path):
    data = np.loadtxt(path, skiprows=1)
    return data[:, 0], data[:, 1]


def main():
    parser = argparse.ArgumentParser(
        description="Plot the final B2 population-weighted UV-Vis envelope "
                    "(neutral + monoanion_b) with matplotlib.")
    parser.add_argument("--neutral-envelope", required=True)
    parser.add_argument("--monoanion-b-envelope", required=True)
    parser.add_argument("--mismatch-rate", type=float, default=0.083)
    parser.add_argument("--mismatch-ci-low", type=float, default=0.0287)
    parser.add_argument("--mismatch-ci-high", type=float, default=0.2183)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    x_n, y_n = load_envelope(args.neutral_envelope)
    x_b, y_b = load_envelope(args.monoanion_b_envelope)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(x_n, y_n, color="#1f77b4", lw=1.8,
            label=(f"Neutral (surrogate + real HSE06, residual "
                    f"state-selection mismatch risk "
                    f"{100*args.mismatch_rate:.1f}% [{100*args.mismatch_ci_low:.1f}%, "
                    f"{100*args.mismatch_ci_high:.1f}%])"))
    ax.plot(x_b, y_b, color="#d62728", lw=1.8,
            label="Monoanion B (100% real HSE06; O···H model valid "
                  "only within a 0.0008 Å training range)")

    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Population-weighted absorption (arb. units)")
    ax.set_title("B2 final UV-Vis envelope, population-weighted")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_xlim(min(x_n.min(), x_b.min()), max(x_n.max(), x_b.max()))
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
