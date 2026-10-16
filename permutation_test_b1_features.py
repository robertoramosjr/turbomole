#!/usr/bin/env python3
"""
permutation_test_b1_features.py

Calibrates the multiple-comparisons risk in select_b1_feature_model.py's
"try 5 features, take the best single-feature Q2" procedure: testing
several features and reporting only the best one's Q2 has a built-in
selection bias -- even pure noise features will sometimes produce a
positive Q2 by chance, and the more features tried, the more likely
one of them looks good for no real reason. This runs the standard
control: shuffle the target (breaking any real feature-target
relationship, keeping the features' own values and correlations
intact), rerun the exact same "best Q2 among N features" procedure on
the shuffled data, and repeat many times to build the null distribution
of "best Q2 you'd see by chance alone, with this many candidate
features and this sample size."

Uses the closed-form (exact, not approximate) leave-one-out formula for
ordinary least squares -- the "PRESS" shortcut via the hat matrix:

    e_loo_i = e_i / (1 - h_ii)

where e_i is the ordinary (full-fit) residual and h_ii is that point's
leverage (diagonal of the hat matrix). This is mathematically identical
to actually refitting n times and predicting the held-out point, just
far cheaper -- needed here because thousands of permutations x several
features x an N-point LOO would otherwise mean fitting the model
tens of thousands of times per state.

Usage:
    python ~/work_turbomole/scripts/permutation_test_b1_features.py \
        --training-data layer_b/b1_training_data.csv \
        --features layer_b/oh_distance_feature.csv \
                   layer_b/planarity_feature.csv \
                   layer_b/population_feature.csv \
        --state-label monoanion_a \
        --n-permutations 1000 \
        --output layer_b/b1_permutation_test.dat
"""

import argparse
import csv

import numpy as np


def load_data(training_path, feature_paths, target_column):
    with open(training_path) as f:
        training = {row["job_dir"]: float(row[target_column])
                    for row in csv.DictReader(f)}

    feature_values = {}
    feature_names = []
    for feat_path in feature_paths:
        with open(feat_path) as f:
            reader = csv.DictReader(f)
            cols = [c for c in reader.fieldnames if c != "job_dir"]
            feature_names.extend(cols)
            for row in reader:
                feature_values.setdefault(row["job_dir"], {})
                for c in cols:
                    feature_values[row["job_dir"]][c] = float(row[c])

    job_dirs = [jd for jd in training if jd in feature_values
                and all(c in feature_values[jd] for c in feature_names)]
    X = np.array([[feature_values[jd][c] for c in feature_names] for jd in job_dirs])
    y = np.array([training[jd] for jd in job_dirs])
    return job_dirs, X, y, feature_names


def loo_q2_closed_form(x, y):
    """Exact leave-one-out Q2 for a single-feature OLS (intercept +
    slope), via the PRESS/hat-matrix shortcut -- no refitting loop."""
    n = len(y)
    design = np.column_stack([np.ones(n), x])
    hat = design @ np.linalg.pinv(design.T @ design) @ design.T
    leverage = np.diag(hat)
    y_pred_full = hat @ y
    resid = y - y_pred_full
    # guard against leverage == 1 (a point that perfectly determines its
    # own fit -- can't happen for n>=3 well-posed simple regression, but
    # guard anyway rather than divide by zero silently)
    denom = np.clip(1 - leverage, 1e-10, None)
    loo_resid = resid / denom
    y_loo_pred = y - loo_resid
    ss_tot = np.sum((y - y.mean()) ** 2)
    ss_res = np.sum(loo_resid ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def best_q2_across_features(X, y):
    return max(loo_q2_closed_form(X[:, i], y) for i in range(X.shape[1]))


def main():
    parser = argparse.ArgumentParser(
        description="Permutation-test calibration for the 'best Q2 among "
                    "several single-feature LOO-OLS models' selection "
                    "procedure used for B1's feature choice.")
    parser.add_argument("--training-data", required=True)
    parser.add_argument("--features", required=True, nargs="+")
    parser.add_argument("--target-column", default="delta_cm1",
                         help="Column name in --training-data holding the "
                              "delta target [default: delta_cm1, B1's "
                              "carbonyl shift -- pass e.g. delta_energy_eV "
                              "or delta_log10_f for B2]")
    parser.add_argument("--state-label", required=True)
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.n_permutations < 500:
        raise SystemExit("--n-permutations must be at least 500.")

    job_dirs, X, y, feature_names = load_data(
        args.training_data, args.features, args.target_column)
    n = len(y)
    print(f"=== {args.state_label} (N={n}, {len(feature_names)} candidate features) ===")

    # sanity check: closed-form LOO should match the refit-based LOO
    # already used elsewhere in this project, for the first feature
    q2_check = loo_q2_closed_form(X[:, 0], y)
    print(f"Sanity check (closed-form LOO Q2 for '{feature_names[0]}'): "
          f"{q2_check:.4f} (compare against select_b1_feature_model.py's "
          f"refit-based value for the same feature -- should match to "
          f"numerical precision)")

    per_feature_q2 = [loo_q2_closed_form(X[:, i], y) for i in range(X.shape[1])]
    observed_best = max(per_feature_q2)
    observed_best_feature = feature_names[int(np.argmax(per_feature_q2))]
    print(f"\nObserved (real target): best Q2 = {observed_best:.4f} "
          f"({observed_best_feature})")
    for name, q in zip(feature_names, per_feature_q2):
        print(f"  {name}: Q2 = {q:.4f}")

    rng = np.random.default_rng(args.seed)
    null_best_q2 = np.zeros(args.n_permutations)
    for k in range(args.n_permutations):
        y_shuffled = rng.permutation(y)
        null_best_q2[k] = best_q2_across_features(X, y_shuffled)

    percentile_observed = float(np.mean(null_best_q2 <= observed_best) * 100)
    p_value = float(np.mean(null_best_q2 >= observed_best))
    q2_95th = float(np.percentile(null_best_q2, 95))

    print(f"\n--- Permutation null distribution ({args.n_permutations} shuffles) ---")
    print(f"  Null 'best-of-{len(feature_names)}' Q2: "
          f"mean={null_best_q2.mean():.4f}, std={null_best_q2.std():.4f}, "
          f"median={np.median(null_best_q2):.4f}")
    print(f"  Null 95th percentile (the 'honest threshold' for future results): "
          f"Q2 = {q2_95th:.4f}")
    print(f"  Observed Q2 = {observed_best:.4f} sits at the "
          f"{percentile_observed:.1f}th percentile of the null distribution")
    print(f"  Empirical p-value (fraction of null >= observed) = {p_value:.4f}")

    if observed_best > 0:
        if percentile_observed >= 95:
            print(f"  => Observed Q2 clears the 95th-percentile null threshold -- "
                  f"NOT easily explained by chance alone across "
                  f"{len(feature_names)} feature tries.")
        elif percentile_observed >= 90:
            print(f"  => Observed Q2 is between the 90th and 95th percentile of "
                  f"the null -- borderline. Flagging explicitly per instructions: "
                  f"do not treat this as safely above the multiple-comparisons "
                  f"noise floor without additional evidence (e.g. more labeled "
                  f"points, or a mechanistic reason to expect this specific "
                  f"feature).")
        else:
            print(f"  => Observed Q2 does NOT clear the 90th percentile of the "
                  f"null -- a result this good is not unusual even with zero "
                  f"real signal, given {len(feature_names)} features were tried "
                  f"and only the best one reported.")
    else:
        # For a non-positive observed result, "typical" means falling in the
        # ordinary body of the null (neither an unusually bad nor, in
        # principle, an unusually good draw) -- only a LOW-tail percentile
        # (worse than nearly all pure-noise draws) is the actual anomaly
        # worth flagging; being above the null's median is not surprising on
        # its own, since a negative-mean null distribution routinely has
        # plenty of draws above its own median.
        if percentile_observed >= 10:
            print(f"  => The 'no signal in any of {len(feature_names)} features' "
                  f"result (percentile {percentile_observed:.1f}) sits within "
                  f"the ordinary body of the null distribution -- consistent "
                  f"with genuinely no learnable signal, not with unusually bad "
                  f"features or a data problem.")
        else:
            print(f"  => NOTE: the observed best-of-{len(feature_names)} Q2 is "
                  f"worse than nearly all pure-noise draws (percentile "
                  f"{percentile_observed:.1f}) -- that is itself unusual and "
                  f"worth checking the features/data for something anomalous, "
                  f"rather than assuming this is an ordinary null result.")

    with open(args.output, "w") as f:
        f.write("descriptor state,n,n_features,n_permutations,observed_best_q2,"
                "observed_best_feature,null_mean,null_std,null_median,"
                "null_95th_percentile,observed_percentile,empirical_pvalue\n")
        f.write(f"{args.state_label} {n} {len(feature_names)} {args.n_permutations} "
                f"{observed_best:.4f} {observed_best_feature} "
                f"{null_best_q2.mean():.4f} {null_best_q2.std():.4f} "
                f"{np.median(null_best_q2):.4f} {q2_95th:.4f} "
                f"{percentile_observed:.2f} {p_value:.4f}\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
