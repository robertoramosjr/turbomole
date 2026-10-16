#!/usr/bin/env python3
"""
apply_b1_model.py

Applies a closed B1 single-feature linear model (delta = slope*feature
+ intercept, e.g. monoanion A's Bader(O)-based model or monoanion B's
O...H-distance-based model, both from
protonation_and_conformer_investigation.md Sec. 9.9/9.12) to a new set
of conformers, flagging every prediction whose feature value falls
OUTSIDE the training range as extrapolation rather than silently
returning a number with the same implied confidence as an interpolated
prediction. A model validated at Q2~0.28-0.59 inside its training range
carries no guarantee whatsoever outside it -- do not trust an
extrapolated value without independent confirmation (e.g. an actual
label for that conformer).

Usage:
    python ~/work_turbomole/scripts/apply_b1_model.py \
        --slope 4.7223 --intercept -73.9249 \
        --valid-min -1.6034 --valid-max -1.1553 \
        --feature-column bader_charge_carbonyl_o \
        --input-features layer_a/some_new_feature.csv \
        --output layer_b/b1_predictions.csv
"""

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(
        description="Apply a closed B1 single-feature linear model with "
                    "explicit extrapolation flagging.")
    parser.add_argument("--slope", type=float, required=True)
    parser.add_argument("--intercept", type=float, required=True)
    parser.add_argument("--valid-min", type=float, required=True,
                         help="Minimum feature value seen in the training set")
    parser.add_argument("--valid-max", type=float, required=True,
                         help="Maximum feature value seen in the training set")
    parser.add_argument("--feature-column", required=True,
                         help="Column name in --input-features to use")
    parser.add_argument("--input-features", required=True,
                         help="CSV with a 'job_dir' column and --feature-column")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.valid_min > args.valid_max:
        raise SystemExit("--valid-min must be <= --valid-max")

    with open(args.input_features) as f:
        rows = list(csv.DictReader(f))

    results = []
    n_extrapolated = 0
    for row in rows:
        jd = row["job_dir"]
        x = float(row[args.feature_column])
        pred = args.slope * x + args.intercept
        in_range = args.valid_min <= x <= args.valid_max
        if not in_range:
            n_extrapolated += 1
        results.append((jd, x, pred, in_range))

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["job_dir", args.feature_column, "predicted_delta_cm1",
                          "in_training_range", "reliability"])
        for jd, x, pred, in_range in results:
            reliability = "interpolation" if in_range else "EXTRAPOLATION -- do not trust without independent confirmation"
            writer.writerow([jd, f"{x:.4f}", f"{pred:.4f}", int(in_range), reliability])

    print(f"Applied model to {len(results)} conformers.")
    print(f"Training range: [{args.valid_min:.4f}, {args.valid_max:.4f}]")
    print(f"{n_extrapolated} / {len(results)} conformers fall OUTSIDE the "
          f"training range (flagged as extrapolation, not to be trusted at "
          f"the same confidence as the Q2 measured inside the range).")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
