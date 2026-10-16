#!/usr/bin/env python3
"""
plot_b2_final_spectrum_veusz.py

Polished Veusz version of the final B2 population-weighted UV-Vis
figure (neutral + monoanion_b), extending plot_optical.py's pattern
(single broadened envelope) to overlay two population-weighted
envelopes with an explicit reliability annotation for each state's
known limitation:
  - neutral: residual state-selection mismatch risk (a categorical
    error, not folded into any error bar) -- see
    protonation_and_conformer_investigation.md Sec. 9.25-9.27.
  - monoanion_b: the B2 O...H model is validated only within a
    0.0008 A training range (Sec. 9.21) -- not an extrapolation risk
    for THIS figure specifically (that state's ensemble is 100% real
    HSE06, no model applied), but noted per the request to reflect it
    in the legend regardless.

NOTE (2026-09-11): veusz is not installed in any of this machine's 7
conda environments -- this script is written to run wherever Veusz IS
available (e.g. the local machine that produced the existing .vsz
files in this project, per their embedded ImportPath). It has not been
executed on this machine. plot_b2_final_spectrum_matplotlib.py is the
fallback that WAS run here, for immediate viewing.

Usage (on a machine with veusz installed):
    python plot_b2_final_spectrum_veusz.py \
        --neutral-envelope b2_final_spectrum_neutral.dat \
        --monoanion-b-envelope b2_final_spectrum_monoanion_b.dat \
        --mismatch-rate 0.083 --mismatch-ci "[2.9%, 21.8%]" \
        --output b2_final_spectrum
"""

import argparse

import veusz.embed as veusz


def main():
    parser = argparse.ArgumentParser(
        description="Plot the final B2 population-weighted UV-Vis envelope "
                    "(neutral + monoanion_b) in Veusz.")
    parser.add_argument("--neutral-envelope", required=True)
    parser.add_argument("--monoanion-b-envelope", required=True)
    parser.add_argument("--mismatch-rate", type=float, default=0.083)
    parser.add_argument("--mismatch-ci", default="[2.9%, 21.8%]")
    parser.add_argument("--output", default="b2_final_spectrum")
    args = parser.parse_args()

    g = veusz.Embedded(args.output, hidden=True)
    g.To(g.Add("page", name="page1"))
    g.Set("width", "16cm")
    g.Set("height", "10cm")
    g.To(g.Add("graph", name="graph1"))

    g.ImportFile(args.neutral_envelope, "", linked=True)
    g.ImportFile(args.monoanion_b_envelope, "", linked=True)

    g.Set("x/label", "Energy (eV)")
    g.Set("y/label", "Population-weighted absorption (arb. units)")

    g.To(g.Add("xy", name="neutral"))
    g.Set("xData", "energy_eV")
    g.Set("yData", "absorption_weighted")
    g.Set("marker", "none")
    g.Set("PlotLine/color", "blue")
    g.Set("key", f"Neutral (residual state-selection mismatch risk "
                  f"{100*args.mismatch_rate:.1f}% {args.mismatch_ci})")
    g.To("..")

    g.To(g.Add("xy", name="monoanion_b"))
    g.Set("xData", "energy_eV_2")  # Veusz auto-suffixes duplicate dataset names on the 2nd ImportFile
    g.Set("yData", "absorption_weighted_2")
    g.Set("marker", "none")
    g.Set("PlotLine/color", "red")
    g.Set("key", "Monoanion B (100% real HSE06; O...H model valid "
                 "only within a 0.0008 A training range)")
    g.To("..")

    g.Add("key", Border__hide=True)

    g.Save(f"{args.output}.vsz")
    g.Export(f"{args.output}.pdf")
    print(f"Wrote {args.output}.vsz and {args.output}.pdf")


if __name__ == "__main__":
    main()
