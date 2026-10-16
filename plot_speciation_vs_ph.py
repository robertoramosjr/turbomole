#!/usr/bin/env python3
"""
plot_speciation_vs_ph.py

Population fraction of each protonation microstate of artepilin C vs.
pH, from microscopic pKa values (protonation_and_conformer_investigation.md
Sec. 9.34). Microstates: N (neutral, COOH/OH), A (carboxylate, COO-/OH),
B (phenolate, COOH/O-), D (dianion). Relative to N, in log10 units:

    w_N = 0
    w_A = pH - pKa_NA
    w_B = pH - pKa_NB
    w_D = 2 pH - pKa_NA - pKa_AD      (route via A)
        = 2 pH - pKa_NB - pKa_BD      (route via B)
          -- if both routes are given, their mean is used and the
             thermodynamic-cycle closure error
             |(pKa_NA+pKa_AD) - (pKa_NB+pKa_BD)| is reported.

fraction_i = 10^w_i / sum_j 10^w_j  (log-sum-exp, numerically stable).

PARTIAL MODE: if pKa_NB / pKa_AD are not given, only N and A are
plotted (a two-state model, valid only where the phenol has not yet
started to ionize -- roughly pH below pKa_phenol - 2) and the figure
says so explicitly. No value is ever invented for a missing pKa.

--band draws the sensitivity to a uniform +/- shift of every input pKa
(NOT a statistical confidence interval -- the calibration has N=1 per
chemical class, so a band size must be chosen by the caller, e.g. the
observed spread between two reference-compound residuals).

Usage:
    python ~/work_turbomole/scripts/plot_speciation_vs_ph.py \
        --pka-na 5.30 --experimental-pka-na 4.65 --band 0.65 \
        --label "carboxyl pKa: calibrated on acetic acid" \
        --output-prefix speciation_partial
    python ~/work_turbomole/scripts/plot_speciation_vs_ph.py \
        --pka-na 5.30 --pka-nb 10.2 --pka-ad 9.9 --output-prefix speciation_full
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt

COLORS = {"N": "#1f77b4", "A": "#2ca02c", "B": "#d62728", "D": "#9467bd"}
LABELS = {"N": "Neutral (COOH, OH)", "A": "Monoanion A (COO$^-$, OH)",
          "B": "Monoanion B (COOH, O$^-$)", "D": "Dianion (COO$^-$, O$^-$)"}


def log_weights(ph, pka_na, pka_nb, pka_ad, pka_bd):
    w = {"N": np.zeros_like(ph), "A": ph - pka_na}
    if pka_nb is not None:
        w["B"] = ph - pka_nb
    if pka_ad is not None or pka_bd is not None:
        routes = []
        if pka_ad is not None:
            routes.append(pka_na + pka_ad)
        if pka_nb is not None and pka_bd is not None:
            routes.append(pka_nb + pka_bd)
        w["D"] = 2 * ph - np.mean(routes)
    return w


def fractions(ph, pka_na, pka_nb, pka_ad, pka_bd):
    w = log_weights(ph, pka_na, pka_nb, pka_ad, pka_bd)
    keys = list(w)
    stack = np.array([w[k] for k in keys])
    stack = stack - stack.max(axis=0)
    p = 10 ** stack
    p = p / p.sum(axis=0)
    return dict(zip(keys, p))


def main():
    parser = argparse.ArgumentParser(
        description="Plot microstate population fractions of artepilin C vs "
                    "pH from microscopic pKa values (partial mode if the "
                    "phenol-related pKa values are not available).")
    parser.add_argument("--pka-na", type=float, required=True,
                         help="N -> A (carboxyl) pKa")
    parser.add_argument("--pka-nb", type=float, default=None,
                         help="N -> B (phenol) pKa")
    parser.add_argument("--pka-ad", type=float, default=None,
                         help="A -> D (phenol, from the carboxylate) pKa")
    parser.add_argument("--pka-bd", type=float, default=None,
                         help="B -> D (carboxyl, from the phenolate) pKa; "
                              "used only for the cycle-closure check/mean")
    parser.add_argument("--experimental-pka-na", type=float, default=None,
                         help="Experimental carboxyl pKa, drawn as a marker "
                              "(Camuri/Pazin: 4.65)")
    parser.add_argument("--band", type=float, default=0.0,
                         help="Uniform +/- pKa shift for a sensitivity band "
                              "(not a confidence interval) [default: 0]")
    parser.add_argument("--label", default="",
                         help="Free text describing where the pKa values came from")
    parser.add_argument("--ph-min", type=float, default=0.0)
    parser.add_argument("--ph-max", type=float, default=14.0)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    ph = np.linspace(args.ph_min, args.ph_max, 1401)
    base = fractions(ph, args.pka_na, args.pka_nb, args.pka_ad, args.pka_bd)
    partial = "B" not in base or "D" not in base

    closure = None
    if None not in (args.pka_nb, args.pka_ad, args.pka_bd):
        closure = abs((args.pka_na + args.pka_ad) - (args.pka_nb + args.pka_bd))
        print(f"Thermodynamic-cycle closure error: {closure:.3f} pKa units "
              f"(route via A vs via B)")

    with open(f"{args.output_prefix}.dat", "w") as f:
        f.write("descriptor pH," + ",".join(f"f_{k}" for k in base) + "\n")
        for i, x in enumerate(ph):
            f.write(f"{x:.3f} " + " ".join(f"{base[k][i]:.6f}" for k in base) + "\n")

    fig, ax = plt.subplots(figsize=(8, 5.2))
    for k, y in base.items():
        ax.plot(ph, y, color=COLORS[k], lw=2, label=LABELS[k])
        if args.band > 0:
            lo = fractions(ph, args.pka_na - args.band,
                           None if args.pka_nb is None else args.pka_nb - args.band,
                           None if args.pka_ad is None else args.pka_ad - args.band,
                           None if args.pka_bd is None else args.pka_bd - args.band)[k]
            hi = fractions(ph, args.pka_na + args.band,
                           None if args.pka_nb is None else args.pka_nb + args.band,
                           None if args.pka_ad is None else args.pka_ad + args.band,
                           None if args.pka_bd is None else args.pka_bd + args.band)[k]
            ax.fill_between(ph, np.minimum(lo, hi), np.maximum(lo, hi),
                            color=COLORS[k], alpha=0.15, lw=0)
    ax.axvline(args.pka_na, color="0.4", ls=":", lw=1)
    ax.text(args.pka_na + 0.12, 0.62, f"pKa(N→A) = {args.pka_na:.2f}",
            rotation=90, va="center", ha="left", fontsize=8, color="0.3")
    if args.experimental_pka_na is not None:
        ax.axvline(args.experimental_pka_na, color="k", ls="--", lw=1)
        ax.text(args.experimental_pka_na - 0.12, 0.62,
                f"exp. pKa = {args.experimental_pka_na:.2f}",
                rotation=90, va="center", ha="right", fontsize=8)
    if "B" not in base:
        ax.axvspan(8.0, args.ph_max, color="0.85", alpha=0.6, lw=0)
        ax.text((8.0 + args.ph_max) / 2, 0.62,
                "phenol pKa not yet calibrated:\ntwo-state model NOT valid here\n"
                "(shading from pH 8 is a heuristic,\nnot a computed value)",
                ha="center", va="center", fontsize=8, color="0.25")
    elif "D" not in base:
        ax.text(0.5, 0.93, "dianion not included (pKa(A→D) missing): "
                "valid only below ~pKa(N→B)", transform=ax.transAxes,
                ha="center", fontsize=8, color="0.25")
    ax.set_xlim(args.ph_min, args.ph_max)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("pH")
    ax.set_ylabel("Population fraction")
    title = "Artepilin C speciation" + (" (PARTIAL: N ⇌ A only)" if partial else "")
    ax.set_title(title, pad=10)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False)
    if args.label:
        fig.text(0.01, 0.01, args.label, fontsize=7, color="0.3", ha="left", va="bottom",
                 wrap=True)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(f"{args.output_prefix}.png", dpi=200)
    print(f"Wrote {args.output_prefix}.png and {args.output_prefix}.dat "
          f"({'partial' if partial else 'full'} model)")


if __name__ == "__main__":
    main()
