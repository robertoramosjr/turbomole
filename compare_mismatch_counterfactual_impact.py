#!/usr/bin/env python3
"""
compare_mismatch_counterfactual_impact.py

Quantifies the maximum plausible impact of an UNDETECTED state-
selection mismatch (the job_0023/job_0035/job_0036 failure signature,
Sec. 9.25-9.26 of protonation_and_conformer_investigation.md) on the
neutral state's population-weighted UV-Vis intensity, by comparing:

  (a) the ACTUAL aggregate intensity (population * f, summed over all
      ranked conformers), using the real, now-confirmed HSE06 values
      for job_0023/job_0035/job_0036 -- what the spectrum actually is;
  (b) a COUNTERFACTUAL aggregate intensity where those same 3
      conformers instead use their pre-confirmation SURROGATE
      prediction (the v2 slope-intercept model applied to their raw
      sTDA max-f value) -- what the spectrum would have looked like
      had these 3 mismatches gone undetected, as ~40%% of the neutral
      state's population still is.

This is not a Gaussian error bar (a state-selection mismatch is a
categorical error -- wrong electronic state, not a noisy magnitude) --
it is reported as a bounded what-if next to the real spectrum, sized
by the ONLY 3 cases where the true impact is actually known.

Usage:
    python ~/work_turbomole/scripts/compare_mismatch_counterfactual_impact.py \
        --population-table layer_a/population_table_dft.csv \
        --predictions-v2 layer_b/b2_osc_predictions_full_ensemble_v2.csv \
        --stda-confidence layer_b/stda_confidence_full_ensemble.csv \
        --known-mismatches job_0023 job_0035 job_0036 \
        --f-model-intercept -0.2027 --f-model-slope 0.5794 \
        --offset-mean-energy 1.2113 \
        --osc-floor 1e-4 \
        --output layer_b/mismatch_counterfactual_impact.dat
"""

import argparse
import csv
import math


def main():
    parser = argparse.ArgumentParser(
        description="Compare actual vs. counterfactual (undetected-mismatch) "
                    "aggregate spectral intensity for the neutral state.")
    parser.add_argument("--population-table", required=True)
    parser.add_argument("--predictions-v2", required=True)
    parser.add_argument("--stda-confidence", required=True,
                         help="characterize_stda_selection_confidence.py "
                              "output, for the pre-confirmation raw sTDA "
                              "max-f value of the known mismatches")
    parser.add_argument("--known-mismatches", nargs="+", required=True)
    parser.add_argument("--f-model-intercept", type=float, required=True)
    parser.add_argument("--f-model-slope", type=float, required=True)
    parser.add_argument("--offset-mean-energy", type=float, required=True)
    parser.add_argument("--osc-floor", type=float, default=1e-4)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.population_table) as f:
        population = {r["job_dir"]: float(r["r2scan3c_population"])
                      for r in csv.DictReader(f)}
    with open(args.predictions_v2) as f:
        pred = {r["job_dir"]: r for r in csv.DictReader(f)}
    with open(args.stda_confidence) as f:
        stda_raw = {r["job_dir"]: (float(r["s1_energy_eV_stda"]), float(r["s1_osc_stda"]))
                    for r in csv.DictReader(f)}

    actual_total = 0.0
    counterfactual_total = 0.0
    per_conformer = []
    for jd, pop in population.items():
        row = pred.get(jd)
        if row is None or row["excluded"] == "True" or not row["predicted_osc"]:
            continue
        f_actual = float(row["predicted_osc"])
        actual_total += pop * f_actual

        if jd in args.known_mismatches:
            e_raw, f_raw = stda_raw[jd]
            log_f_raw = math.log10(max(f_raw, args.osc_floor))
            f_counterfactual = 10 ** (args.f_model_intercept + args.f_model_slope * log_f_raw)
            counterfactual_total += pop * f_counterfactual
            per_conformer.append((jd, pop, f_actual, f_counterfactual))
        else:
            counterfactual_total += pop * f_actual

    print(f"=== Known mismatches: per-conformer actual vs. counterfactual ===")
    for jd, pop, f_actual, f_cf in per_conformer:
        pct_pop = 100 * pop
        print(f"  {jd}: population={pop:.4f} ({pct_pop:.2f}%), "
              f"f_actual(real)={f_actual:.4f}, "
              f"f_counterfactual(undetected surrogate)={f_cf:.4f}, "
              f"ratio={f_cf / f_actual:.1f}x" if f_actual > 1e-9 else "")

    delta = counterfactual_total - actual_total
    pct_shift = 100 * delta / actual_total if actual_total > 0 else float("nan")
    print(f"\n=== Aggregate population-weighted intensity (sum of pop x f) ===")
    print(f"  Actual (real values for the 3 known mismatches):        {actual_total:.5f}")
    print(f"  Counterfactual (undetected, surrogate for the 3):       {counterfactual_total:.5f}")
    print(f"  Difference: {delta:+.5f} ({pct_shift:+.2f}% of the actual aggregate)")
    print(f"\n  This bounds the plausible aggregate distortion from 3 "
          f"undetected mismatches of this specific weight/magnitude "
          f"profile -- NOT a general error bar for the remaining ~8.3% "
          f"background rate among still-unconfirmed conformers, which "
          f"could in principle be larger or smaller depending on which "
          f"conformers are affected.")

    with open(args.output, "w") as f:
        f.write("descriptor job_dir,population,f_actual,f_counterfactual\n")
        for jd, pop, f_actual, f_cf in per_conformer:
            f.write(f"{jd} {pop:.4f} {f_actual:.4f} {f_cf:.4f}\n")
        f.write(f"\ndescriptor aggregate_actual,aggregate_counterfactual,"
                f"delta,pct_shift\n")
        f.write(f"{actual_total:.5f} {counterfactual_total:.5f} {delta:.5f} {pct_shift:.2f}\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
