#!/usr/bin/env python3
"""
build_b2_osc_delta_target.py

Builds B2's delta-learning target table by joining the HSE06 S1
reference (batch_extract_s1_hse06.py) with the character-matched sTDA
proxy (batch_stda_s1_character_matched.py): delta_energy_eV = HSE06 -
sTDA (energy), delta_log10_f = log10(HSE06 f) - log10(sTDA f), both in
the flat units select_b1_feature_model.py / permutation_test_b1_features.py
expect (one row per job_dir, plus a --target-column-selectable value).

Oscillator strengths are floored at --osc-floor before the log10 (both
sides), matching the "near-zero" convention already used in
diagnose_stda_oscillator.py (1e-3) but defaulting slightly lower
(1e-4) since several genuinely dark states in this project's own data
sit right at 1e-4 (see monoanion_a, dianion) -- flooring further below
that would treat real near-dark labels as if they were noise-free
exact zeros, and flooring above it would clip real signal.

Usage:
    python ~/work_turbomole/scripts/build_b2_osc_delta_target.py \
        --hse06 layer_b/s1_hse06.csv \
        --stda layer_b/s1_stda_character_matched.csv \
        --osc-floor 1e-4 \
        --output layer_b/b2_osc_delta_target.csv
"""

import argparse
import csv
import math


def main():
    parser = argparse.ArgumentParser(
        description="Build B2's energy + log(f) delta target table from "
                    "HSE06 and character-matched sTDA S1 CSVs.")
    parser.add_argument("--hse06", required=True)
    parser.add_argument("--stda", required=True)
    parser.add_argument("--osc-floor", type=float, default=1e-4)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.hse06) as f:
        hse06 = {row["job_dir"]: (float(row["s1_energy_eV"]), float(row["s1_osc_len"]))
                  for row in csv.DictReader(f)}
    with open(args.stda) as f:
        stda = {row["job_dir"]: (float(row["s1_energy_eV_stda"]), float(row["s1_osc_stda"]))
                for row in csv.DictReader(f)}

    common = sorted(set(hse06) & set(stda))
    rows = []
    for jd in common:
        e_h, f_h = hse06[jd]
        e_s, f_s = stda[jd]
        f_h_floored = max(f_h, args.osc_floor)
        f_s_floored = max(f_s, args.osc_floor)
        rows.append({
            "job_dir": jd,
            "delta_energy_eV": e_h - e_s,
            "delta_log10_f": math.log10(f_h_floored) - math.log10(f_s_floored),
            "hse06_energy_eV": e_h, "hse06_osc": f_h,
            "stda_energy_eV": e_s, "stda_osc": f_s,
        })

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Built {len(rows)} / {len(set(hse06) | set(stda))} joined rows.")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
