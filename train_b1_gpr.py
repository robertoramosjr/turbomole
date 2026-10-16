#!/usr/bin/env python3
"""
train_b1_gpr.py

Trains the B1 delta-learning surrogate (GPR, Matern-5/2 kernel,
predicting HSE06-minus-GFN2xTB correction to the C=O stretch
frequency) on one protonation state's labeled subset, and evaluates it
by leave-one-out cross-validation (appropriate here: only 11-20
labeled points per state, a single held-out split would be too noisy
to trust). This is the active-learning loop's "iteration 0" bootstrap
-- not yet actually active learning (no new labels are picked here),
just the starting point whose ranking metric the loop will try to
improve.

Feature: O(carbonyl)...H(nearest donor) distance (compute_b1_features.py)
-- the one Paper 1 flagged as correlating with the band shift. Reported
here whether it actually carries useful signal for THIS labeled sample
before assuming it does (per this project's documented caution: "vale
verificar se essa unica feature geometrica ja explica a maior parte da
variancia antes de depender so do GPR sobre diedros").

Ranking metric: Spearman rho between true and leave-one-out-predicted
delta, over the labeled subset. Recall@10 is not computed at this
stage -- with only 11-20 total labeled conformers, "top 10" is most or
all of the sample and the metric is not meaningful; it becomes
meaningful once the surrogate ranks the full ~150-340-conformer Layer A
ensemble, which is the next stage, not this one.

Usage:
    python ~/work_turbomole/scripts/train_b1_gpr.py \
        --training-data layer_b/b1_training_data.csv \
        --features layer_b/oh_distance_feature.csv \
        --output layer_b/b1_loo_predictions.csv
"""

import argparse
import csv

import numpy as np
from scipy.stats import spearmanr
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler


def main():
    parser = argparse.ArgumentParser(
        description="Train + leave-one-out-evaluate the B1 GPR delta-learning "
                    "surrogate for one protonation state.")
    parser.add_argument("--training-data", required=True,
                         help="batch_hse06_co_stretch.py output "
                              "(b1_training_data.csv)")
    parser.add_argument("--features", required=True, nargs="+",
                         help="One or more feature CSVs, each with a "
                              "'job_dir' column plus one or more numeric "
                              "feature columns (e.g. oh_distance_feature.csv, "
                              "dihedral_feature.csv). All non-job_dir columns "
                              "across all files are used as GPR input features.")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.training_data) as f:
        training = {row["job_dir"]: float(row["delta_cm1"])
                    for row in csv.DictReader(f)}

    feature_values = {}   # job_dir -> {col_name: value}
    feature_names = []
    for feat_path in args.features:
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
    missing = set(training) - set(job_dirs)
    if missing:
        print(f"WARNING: {len(missing)} labeled conformers have incomplete "
              f"feature values, excluded: {sorted(missing)}")

    X = np.array([[feature_values[jd][c] for c in feature_names] for jd in job_dirs])
    y = np.array([training[jd] for jd in job_dirs])
    n = len(job_dirs)

    print(f"{n} labeled conformers. Features: {feature_names}")
    for i, c in enumerate(feature_names):
        print(f"  {c}: range {X[:, i].min():.3f} - {X[:, i].max():.3f}")
    print(f"Target (delta) range: {y.min():.2f} to {y.max():.2f} cm-1.")

    if n < 4:
        print("Too few points for leave-one-out cross-validation "
              "(need at least a handful) -- reporting the raw target "
              "statistics only, no surrogate evaluation.")
        print(f"Mean delta: {y.mean():.2f} cm-1, std: {y.std():.2f} cm-1")
        with open(args.output, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["job_dir", "true_delta_cm1", "loo_predicted_delta_cm1"])
        return

    x_scaler = StandardScaler()
    y_scaler = StandardScaler()
    X_scaled_full = x_scaler.fit_transform(X)
    y_scaled_full = y_scaler.fit_transform(y.reshape(-1, 1)).ravel()

    kernel = Matern(length_scale=1.0, nu=2.5) + WhiteKernel(noise_level=0.1)

    loo = LeaveOneOut()
    predictions = np.zeros(n)
    std_devs = np.zeros(n)
    for train_idx, test_idx in loo.split(X):
        x_scaler_fold = StandardScaler().fit(X[train_idx])
        y_scaler_fold = StandardScaler().fit(y[train_idx].reshape(-1, 1))
        X_train = x_scaler_fold.transform(X[train_idx])
        y_train = y_scaler_fold.transform(y[train_idx].reshape(-1, 1)).ravel()
        X_test = x_scaler_fold.transform(X[test_idx])

        gpr = GaussianProcessRegressor(kernel=kernel, normalize_y=False,
                                        n_restarts_optimizer=5, random_state=0)
        gpr.fit(X_train, y_train)
        pred_scaled, std_scaled = gpr.predict(X_test, return_std=True)
        predictions[test_idx] = y_scaler_fold.inverse_transform(
            pred_scaled.reshape(-1, 1)).ravel()
        std_devs[test_idx] = std_scaled * y_scaler_fold.scale_

    rho, pval = spearmanr(y, predictions)
    mae = np.mean(np.abs(y - predictions))
    baseline_mae = np.mean(np.abs(y - y.mean()))  # predicting the mean for everyone

    print(f"\nLeave-one-out results:")
    print(f"  Spearman rho (true vs. predicted delta): {rho:.3f} (p={pval:.3f})")
    print(f"  MAE: {mae:.2f} cm-1 (vs. {baseline_mae:.2f} cm-1 for predicting "
          f"the per-state mean delta for every conformer)")
    if mae >= baseline_mae:
        print("  NOTE: the GPR does not beat just predicting the mean delta -- "
              "the O...H distance feature alone may not carry enough signal "
              "for this labeled sample. Worth checking a raw feature-vs-delta "
              "scatter before adding more features.")

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["job_dir", *feature_names, "true_delta_cm1",
                          "loo_predicted_delta_cm1", "loo_predicted_std_cm1"])
        for jd, xrow, yt, yp, ys in zip(job_dirs, X, y, predictions, std_devs):
            writer.writerow([jd, *[f"{v:.4f}" for v in xrow],
                              f"{yt:.3f}", f"{yp:.3f}", f"{ys:.3f}"])
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
