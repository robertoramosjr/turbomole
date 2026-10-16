#!/usr/bin/env python3
"""
diagnose_stda_oscillator.py

Diagnoses why the HSE06/sTDA oscillator-strength RATIO is wildly
uncontrolled (1.5x-235x for neutral, 140x-8000x for monoanion B, per
batch_stda_s1.py's initial ax=0.0 run) despite the S1 excitation
ENERGY delta being well-behaved. A ratio blows up whenever the
denominator (f_sTDA) is small, regardless of whether the two methods
actually disagree by much in absolute terms -- this script checks
whether that is exactly what is happening (the extreme ratios coming
from near-zero f_sTDA transitions, where even a tiny absolute
difference produces a huge ratio) rather than a genuinely uncontrolled
absolute discrepancy.

Reports, per state, into a .dat table:
  - log10(f_HSE06) vs log10(f_sTDA) pairs (for an external log-log
    scatter plot)
  - the ABSOLUTE error |f_HSE06 - f_sTDA| per conformer (not the ratio)
  - how many conformers have f_sTDA < 1e-3, and whether the most
    extreme ratios (>100x) are concentrated in that subset

Usage:
    python ~/work_turbomole/scripts/diagnose_stda_oscillator.py \
        --hse06 layer_b/s1_hse06.csv --stda layer_b/s1_stda.csv \
        --state-label neutral --output layer_b/oscillator_diagnostic.dat
"""

import argparse
import csv
import math


def main():
    parser = argparse.ArgumentParser(
        description="Diagnose HSE06 vs. sTDA oscillator strength "
                    "discrepancy: absolute error vs. ratio, and whether "
                    "extreme ratios trace back to near-zero f_sTDA.")
    parser.add_argument("--hse06", required=True, help="batch_extract_s1_hse06.py output")
    parser.add_argument("--stda", required=True, help="batch_stda_s1.py output")
    parser.add_argument("--state-label", required=True)
    parser.add_argument("--near-zero-threshold", type=float, default=1e-3,
                         help="f_sTDA below this counts as 'near zero' [default: 1e-3]")
    parser.add_argument("--extreme-ratio-threshold", type=float, default=100.0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.hse06) as f:
        hse06 = {row["job_dir"]: float(row["s1_osc_len"]) for row in csv.DictReader(f)}
    with open(args.stda) as f:
        stda = {row["job_dir"]: float(row["s1_osc_stda"]) for row in csv.DictReader(f)}

    common = sorted(set(hse06) & set(stda))
    n = len(common)
    print(f"=== {args.state_label} (N={n}) ===")

    rows = []
    n_near_zero = 0
    n_extreme_ratio = 0
    n_extreme_and_near_zero = 0
    abs_errors = []
    for jd in common:
        f_h, f_s = hse06[jd], stda[jd]
        abs_err = abs(f_h - f_s)
        ratio = f_h / f_s if f_s > 1e-12 else float("inf")
        is_near_zero = f_s < args.near_zero_threshold
        is_extreme = ratio > args.extreme_ratio_threshold
        if is_near_zero:
            n_near_zero += 1
        if is_extreme:
            n_extreme_ratio += 1
            if is_near_zero:
                n_extreme_and_near_zero += 1
        abs_errors.append(abs_err)
        log_h = math.log10(f_h) if f_h > 1e-12 else float("-inf")
        log_s = math.log10(f_s) if f_s > 1e-12 else float("-inf")
        rows.append((jd, f_h, f_s, abs_err, ratio, log_h, log_s, is_near_zero, is_extreme))

    print(f"\nAbsolute error |f_HSE06 - f_sTDA|: "
          f"mean={sum(abs_errors)/n:.4f}, max={max(abs_errors):.4f}, "
          f"min={min(abs_errors):.4f}")
    print(f"Conformers with f_sTDA < {args.near_zero_threshold}: "
          f"{n_near_zero} / {n}")
    print(f"Conformers with ratio > {args.extreme_ratio_threshold}x: "
          f"{n_extreme_ratio} / {n}")
    if n_extreme_ratio > 0:
        pct = 100 * n_extreme_and_near_zero / n_extreme_ratio
        print(f"Of those extreme-ratio conformers, {n_extreme_and_near_zero} "
              f"({pct:.0f}%) also have near-zero f_sTDA")
        if pct >= 80:
            print("=> Extreme ratios are concentrated in the near-zero-f_sTDA "
                  "subset -- consistent with a denominator artifact (tiny "
                  "denominator inflating an otherwise modest absolute "
                  "mismatch), not necessarily an uncontrolled absolute "
                  "discrepancy across the board.")
        else:
            print("=> Extreme ratios are NOT mostly explained by near-zero "
                  "f_sTDA -- some conformers disagree by orders of magnitude "
                  "in ABSOLUTE terms too, which the near-zero-denominator "
                  "explanation does not cover.")
    else:
        print("No extreme-ratio conformers at this threshold.")

    with open(args.output, "w") as f:
        f.write("descriptor job_dir,f_hse06,f_stda,abs_error,ratio,log10_f_hse06,"
                "log10_f_stda,is_near_zero,is_extreme_ratio\n")
        for jd, f_h, f_s, ae, r, lh, ls, nz, ex in rows:
            r_str = f"{r:.4f}" if r != float("inf") else "inf"
            lh_str = f"{lh:.4f}" if lh != float("-inf") else "-inf"
            ls_str = f"{ls:.4f}" if ls != float("-inf") else "-inf"
            f.write(f"{jd} {f_h:.6f} {f_s:.6f} {ae:.6f} {r_str} {lh_str} {ls_str} "
                    f"{int(nz)} {int(ex)}\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
