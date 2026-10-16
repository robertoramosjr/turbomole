#!/usr/bin/env python3
"""
test_b1_linear_baseline.py

Tests whether the GPR's negative Q2 on monoanion B (N=11) despite a
significant single-feature Pearson correlation (r=0.839 for the O...H
distance, per audit_b1_features.py) is caused by GPR overparameterization
(4 hyperparameters -- 2 length scales, signal variance, noise level --
fit on 11 points) rather than a genuine absence of learnable signal.

Three checks, all on the SAME leave-one-out (LOO) protocol used for
the GPR (diagnose_b1_surrogate.py / train_b1_gpr.py), so results are
directly comparable:

  1. Simple OLS (1 feature, 2 parameters: slope + intercept) under LOO
     -- if THIS beats the GPR's Q2 with an order of magnitude fewer
     free parameters, overparameterization is the likely explanation.

  2. Leverage/robustness check on the raw (non-LOO) Pearson correlation
     itself: recomputes r after removing each point one at a time (a
     "delete-one" sensitivity analysis of r, distinct from the model's
     own LOO predictions) and Cook's distance for every point in the
     full-data OLS fit, to catch a correlation propped up by a single
     influential point rather than a real trend across the sample.

  3. The same LOO-OLS test repeated for the second feature alone, and
     for both features together (linear, no interaction term) -- to
     see whether a positive Q2 (if any) survives feature choice or is
     an accident of picking exactly the O...H distance.

Usage:
    python ~/work_turbomole/scripts/test_b1_linear_baseline.py \
        --training-data layer_b/b1_training_data.csv \
        --features layer_b/oh_distance_feature.csv layer_b/planarity_feature.csv \
        --state-label monoanion_b \
        --output layer_b/b1_linear_baseline.dat
"""

import argparse
import csv
import itertools

import numpy as np
from scipy.stats import pearsonr


def load_data(training_path, feature_paths):
    with open(training_path) as f:
        training = {row["job_dir"]: float(row["delta_cm1"])
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


def ols_fit_predict(X_train, y_train, x_test):
    """Ordinary least squares with an intercept; X_train is (n, p)."""
    design = np.column_stack([np.ones(len(X_train)), X_train])
    coef, *_ = np.linalg.lstsq(design, y_train, rcond=None)
    x_test_design = np.concatenate([[1.0], x_test])
    return float(x_test_design @ coef)


def loo_ols(X, y):
    n = len(y)
    preds = np.zeros(n)
    for i in range(n):
        train_idx = [j for j in range(n) if j != i]
        preds[i] = ols_fit_predict(X[train_idx], y[train_idx], X[i])
    return preds


def q2_and_mae(y, preds):
    ss_tot = np.sum((y - y.mean()) ** 2)
    ss_res = np.sum((y - preds) ** 2)
    q2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    mae = np.mean(np.abs(y - preds))
    return q2, mae


def cooks_distance(X, y):
    """Cook's distance for every point in a full-data OLS fit
    (intercept + len(feature) slopes)."""
    n = len(y)
    design = np.column_stack([np.ones(n), X])
    p = design.shape[1]
    hat = design @ np.linalg.pinv(design.T @ design) @ design.T
    leverage = np.diag(hat)
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    residuals = y - design @ coef
    mse = np.sum(residuals ** 2) / (n - p)
    with np.errstate(divide="ignore", invalid="ignore"):
        cooks_d = (residuals ** 2 / (p * mse)) * (leverage / (1 - leverage) ** 2)
    return cooks_d, leverage


def delete_one_correlation(x, y):
    """Recompute Pearson r after removing each point one at a time."""
    n = len(y)
    rs = np.zeros(n)
    for i in range(n):
        keep = [j for j in range(n) if j != i]
        r, _ = pearsonr(x[keep], y[keep])
        rs[i] = r
    return rs


def report_single_feature_test(name, x, y, out_lines):
    r_full, p_full = pearsonr(x, y)
    preds = loo_ols(x.reshape(-1, 1), y)
    q2, mae = q2_and_mae(y, preds)
    cooks_d, leverage = cooks_distance(x.reshape(-1, 1), y)
    delete_one_r = delete_one_correlation(x, y)

    print(f"\n--- Single feature: {name} ---")
    print(f"  Full-sample Pearson r = {r_full:.4f} (p={p_full:.4f})")
    print(f"  LOO-OLS (2 params: slope+intercept): Q2 = {q2:.4f}, MAE = {mae:.3f} cm-1")
    print(f"  Delete-one r range: [{delete_one_r.min():.4f}, {delete_one_r.max():.4f}], "
          f"std = {delete_one_r.std():.4f}")
    max_cook_idx = int(np.argmax(cooks_d))
    print(f"  Cook's distance: max = {cooks_d[max_cook_idx]:.4f} at point index {max_cook_idx} "
          f"(rule of thumb: > 4/n = {4 / len(y):.4f} is worth a look)")
    influential = np.where(cooks_d > 4 / len(y))[0]
    if len(influential) > 0:
        print(f"  WARNING: point(s) {list(influential)} exceed the 4/n Cook's "
              f"distance threshold -- check whether removing them changes the "
              f"conclusion (see delete-one r range above for how much r moves).")
    else:
        print("  No single point exceeds the 4/n Cook's distance rule of thumb.")

    out_lines.append(f"{name} {r_full:.4f} {p_full:.4f} {q2:.4f} {mae:.4f} "
                      f"{delete_one_r.min():.4f} {delete_one_r.max():.4f} "
                      f"{delete_one_r.std():.4f} {cooks_d.max():.4f} {max_cook_idx}")
    return q2


def main():
    parser = argparse.ArgumentParser(
        description="Test a simple OLS baseline (far fewer parameters than "
                    "the GPR) plus leverage/Cook's-distance diagnostics, to "
                    "check whether a negative GPR Q2 despite a significant "
                    "single-feature correlation is overparameterization.")
    parser.add_argument("--training-data", required=True)
    parser.add_argument("--features", required=True, nargs="+")
    parser.add_argument("--state-label", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    job_dirs, X, y, feature_names = load_data(args.training_data, args.features)
    n = len(y)
    print(f"=== {args.state_label} (N={n}) ===")

    if n < 4:
        print("N < 4: too few points for this analysis -- skipping "
              "(this script is not meant for the dianion, N=3).")
        with open(args.output, "w") as f:
            f.write(f"# {args.state_label}: N={n}, skipped\n")
        return

    out_lines = []
    q2_per_feature = {}
    for i, name in enumerate(feature_names):
        q2_per_feature[name] = report_single_feature_test(name, X[:, i], y, out_lines)

    print(f"\n--- Combined (all {len(feature_names)} features, linear, no interaction) ---")
    preds_combined = loo_ols(X, y)
    q2_combined, mae_combined = q2_and_mae(y, preds_combined)
    print(f"  LOO-OLS ({1 + len(feature_names)} params): "
          f"Q2 = {q2_combined:.4f}, MAE = {mae_combined:.3f} cm-1")
    cooks_d_combined, _ = cooks_distance(X, y)
    out_lines.append(f"combined_linear NA NA {q2_combined:.4f} {mae_combined:.4f} "
                      f"NA NA NA {cooks_d_combined.max():.4f} "
                      f"{int(np.argmax(cooks_d_combined))}")

    print(f"\n--- Summary: best single-feature Q2 vs. combined ---")
    best_feature = max(q2_per_feature, key=q2_per_feature.get)
    print(f"  Best single feature: {best_feature} (Q2 = {q2_per_feature[best_feature]:.4f})")
    print(f"  Combined linear model: Q2 = {q2_combined:.4f}")
    print(f"  (compare both to the GPR's Q2 from diagnose_b1_surrogate.py for "
          f"this same state -- if a 2-parameter OLS beats a 4+-parameter GPR "
          f"here, overparameterization is the likely explanation for the GPR's "
          f"negative Q2, not absence of signal)")

    with open(args.output, "w") as f:
        f.write("descriptor feature,pearson_r,pearson_p,q2_loo_ols,mae_loo_ols,"
                "delete_one_r_min,delete_one_r_max,delete_one_r_std,"
                "cooks_d_max,cooks_d_max_index\n")
        for line in out_lines:
            f.write(line + "\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
