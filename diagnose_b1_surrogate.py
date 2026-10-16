#!/usr/bin/env python3
"""
diagnose_b1_surrogate.py

Diagnostic for the B1 delta-learning GPR surrogate, run BEFORE trusting
any leave-one-out (LOO) Spearman rho as evidence of real ranking skill.

Motivation: a predictor that always returns the mean of the OTHER
training points is not a null model with rho near 0 under LOO -- it is
an exact linear (decreasing) function of the true value:

    pred_i = mean(y_train_excluding_i) = (sum(y) - y_i) / (n - 1)
           = [n/(n-1)] * mean(y) - [1/(n-1)] * y_i

Since pred_i is a strictly decreasing linear function of y_i, Spearman
rho(y, pred) = -1 EXACTLY (ties aside), for ANY dataset, regardless of
whether the model learned anything. A GPR whose kernel signal is
swamped by noise (or whose features carry no real signal) will fall
back toward predicting close to the training mean in each fold, and
can therefore show a strongly negative LOO rho that looks like "the
model works (just inverted)" but is actually this same mean-reversion
artifact. This script checks for exactly that before task 2 (a joint
multi-state model) is attempted.

Reports, per protonation state, into a Veusz-style .dat table:
  - N, Var(y), std(y), range(y)
  - Var(y_pred_LOO) / Var(y)          -- near 0 means predictions
    collapsed to (near-)constant, a red flag on its own
  - Q² = 1 - SS_res_LOO / SS_tot      -- skill vs. the mean baseline;
    Q² <= 0 means the model is no better (or worse) than predicting
    the mean, no matter what rho says
  - MAE_LOO and MAE of the "always predict mean-of-others" baseline
  - fitted hyperparameters (length scale per feature -- an ARD Matern
    kernel is used so each feature gets its own value --, signal
    variance, noise level, and the noise/signal ratio) from a SINGLE
    fit on the full labeled set (representative values; actual LOO
    evaluation below refits independently in every fold, which is
    confirmed explicitly in the printed output, not assumed)
  - the mean-of-others control's own Spearman rho (should print at or
    very near -1.000 for every state -- if the real model's rho is
    close to this control's rho, that is evidence for the artifact,
    not for real feature-driven ranking skill)

Usage:
    python ~/work_turbomole/scripts/diagnose_b1_surrogate.py \
        --training-data layer_b/b1_training_data.csv \
        --features layer_b/oh_distance_feature.csv layer_b/planarity_feature.csv \
        --state-label neutral \
        --output layer_b/b1_diagnostic.dat
"""

import argparse
import csv

import numpy as np
from scipy.stats import spearmanr
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.preprocessing import StandardScaler


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


def make_kernel(n_features):
    # ARD (one length scale per feature, not a single shared scalar) so
    # each feature's relevance can be inspected individually; an
    # explicit ConstantKernel makes signal variance its own learnable
    # parameter (sklearn's Matern alone is unit-variance) instead of
    # conflating it with the length scales; a separate WhiteKernel term
    # models observation noise.
    return (ConstantKernel(1.0, (1e-3, 1e3))
            * Matern(length_scale=np.ones(n_features),
                      length_scale_bounds=(1e-3, 1e3), nu=2.5)
            + WhiteKernel(noise_level=1.0, noise_level_bounds=(1e-5, 1e2)))


def loo_predict(X, y, kernel_factory):
    """Leave-one-out predictions, refitting hyperparameters independently
    in every fold (confirmed here, not assumed) -- returns predictions
    array and, per fold, the fitted length scales / noise level / signal
    variance for later inspection of fold-to-fold stability."""
    n = len(y)
    preds = np.zeros(n)
    fold_hparams = []
    for i in range(n):
        train_idx = [j for j in range(n) if j != i]
        x_scaler = StandardScaler().fit(X[train_idx])
        y_scaler = StandardScaler().fit(y[train_idx].reshape(-1, 1))
        X_train = x_scaler.transform(X[train_idx])
        y_train = y_scaler.transform(y[train_idx].reshape(-1, 1)).ravel()
        X_test = x_scaler.transform(X[i:i + 1])

        gpr = GaussianProcessRegressor(kernel=kernel_factory(X.shape[1]),
                                        normalize_y=False,
                                        n_restarts_optimizer=5, random_state=0)
        gpr.fit(X_train, y_train)
        pred_scaled = gpr.predict(X_test)
        preds[i] = y_scaler.inverse_transform(pred_scaled.reshape(-1, 1)).ravel()[0]

        k = gpr.kernel_
        const_part, matern_part = k.k1.k1, k.k1.k2
        white_part = k.k2
        fold_hparams.append({
            "length_scale": np.atleast_1d(matern_part.length_scale).copy(),
            "signal_variance": const_part.constant_value,
            "noise_level": white_part.noise_level,
        })
    return preds, fold_hparams


def mean_of_others_control(y):
    """The degenerate 'model': predict the mean of every OTHER training
    point. Exercises the mean-reversion artifact directly, with no
    features and no fitting at all."""
    n = len(y)
    total = y.sum()
    return np.array([(total - y[i]) / (n - 1) for i in range(n)])


def main():
    parser = argparse.ArgumentParser(
        description="Diagnose whether a B1 GPR's LOO Spearman rho reflects "
                    "real ranking skill or the mean-reversion LOO artifact.")
    parser.add_argument("--training-data", required=True)
    parser.add_argument("--features", required=True, nargs="+")
    parser.add_argument("--state-label", required=True,
                         help="Label for this protonation state (for the report)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    job_dirs, X, y, feature_names = load_data(args.training_data, args.features)
    n = len(y)

    print(f"=== {args.state_label} (N={n}) ===")
    print(f"Features: {feature_names}")

    if n < 3:
        print("N < 3: not enough points for any of this -- skipping.")
        with open(args.output, "w") as f:
            f.write(f"# {args.state_label}: N={n}, insufficient for diagnostics\n")
        return

    ss_tot = np.sum((y - y.mean()) ** 2)
    var_y = np.var(y, ddof=1) if n > 1 else 0.0
    print(f"y (delta_cm1): mean={y.mean():.3f}, std={np.sqrt(var_y):.3f}, "
          f"range=[{y.min():.3f}, {y.max():.3f}] ({y.max() - y.min():.3f} cm-1 wide)")

    # control: mean-of-others, no features, no fitting
    control_preds = mean_of_others_control(y)
    control_rho, control_p = (spearmanr(y, control_preds) if n >= 3 else (np.nan, np.nan))
    control_mae = np.mean(np.abs(y - control_preds))
    print(f"\nControl (predict mean-of-others, no model at all):")
    print(f"  Spearman rho = {control_rho:.4f} (expected: exactly or very near -1.000 "
          f"by construction, ties aside)")
    print(f"  MAE = {control_mae:.3f} cm-1")

    # actual GPR, LOO, hyperparameters refit independently every fold
    # (this loop does exactly that -- see loo_predict's docstring)
    preds, fold_hparams = loo_predict(X, y, make_kernel)
    rho, pval = spearmanr(y, preds)
    ss_res = np.sum((y - preds) ** 2)
    q2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    mae_loo = np.mean(np.abs(y - preds))
    var_ratio = np.var(preds, ddof=1) / var_y if var_y > 0 else float("nan")

    print(f"\nGPR (ARD Matern 5/2 + WhiteKernel), LOO "
          f"(hyperparameters refit independently in every fold):")
    print(f"  Spearman rho = {rho:.4f} (p={pval:.4f})")
    print(f"  Q2 (skill vs. mean baseline) = {q2:.4f} "
          f"({'WORSE than baseline' if q2 < 0 else 'better than baseline' if q2 > 0 else 'ties baseline'})")
    print(f"  MAE_LOO = {mae_loo:.3f} cm-1 (control MAE = {control_mae:.3f} cm-1)")
    print(f"  Var(pred_LOO)/Var(y) = {var_ratio:.4f} "
          f"({'predictions nearly constant -- collapsed toward the mean' if var_ratio < 0.1 else 'predictions retain real spread'})")

    print(f"\n  Rho comparison: model rho = {rho:.4f} vs. mean-reversion "
          f"control rho = {control_rho:.4f}")
    if abs(rho - control_rho) < 0.1:
        print("  => Model rho is essentially indistinguishable from the "
              "mean-reversion artifact -- DO NOT treat this rho as evidence "
              "of real ranking skill.")
    else:
        print("  => Model rho differs meaningfully from the mean-reversion "
              "control -- some real feature-driven signal beyond pure "
              "mean-reversion, though Q2/variance-ratio above still matter "
              "for how much.")

    # representative hyperparameters: one fit on the full labeled set
    x_scaler = StandardScaler().fit(X)
    y_scaler = StandardScaler().fit(y.reshape(-1, 1))
    gpr_full = GaussianProcessRegressor(kernel=make_kernel(X.shape[1]),
                                         normalize_y=False,
                                         n_restarts_optimizer=5, random_state=0)
    gpr_full.fit(x_scaler.transform(X), y_scaler.transform(y.reshape(-1, 1)).ravel())
    k_full = gpr_full.kernel_
    signal_variance = k_full.k1.k1.constant_value
    length_scales = np.atleast_1d(k_full.k1.k2.length_scale)
    noise_level = k_full.k2.noise_level
    noise_signal_ratio = noise_level / signal_variance if signal_variance > 0 else float("inf")
    print(f"\nRepresentative hyperparameters (single fit on the full "
          f"labeled set, standardized features/target; per-fold LOO fits "
          f"above are independent of this):")
    for name, ls in zip(feature_names, length_scales):
        print(f"  length_scale[{name}] = {ls:.4f}")
    print(f"  signal_variance (ConstantKernel) = {signal_variance:.4f}")
    print(f"  noise_level (WhiteKernel) = {noise_level:.4f}")
    print(f"  noise/signal ratio = {noise_signal_ratio:.4f} "
          f"({'noise-dominated -- little room for the kernel to explain '
             'anything beyond a constant' if noise_signal_ratio > 1 else 'signal-dominated'})")

    fold_ls_spread = np.array([fh["length_scale"] for fh in fold_hparams])
    fold_noise_spread = np.array([fh["noise_level"] for fh in fold_hparams])
    print(f"\nFold-to-fold hyperparameter stability (min/max across the "
          f"{n} independent LOO refits):")
    for i, name in enumerate(feature_names):
        print(f"  length_scale[{name}]: [{fold_ls_spread[:, i].min():.4f}, "
              f"{fold_ls_spread[:, i].max():.4f}]")
    print(f"  noise_level: [{fold_noise_spread.min():.4f}, {fold_noise_spread.max():.4f}]")

    with open(args.output, "w") as f:
        f.write("descriptor state,N,var_y,std_y,range_y,control_rho,control_mae,"
                "model_rho,model_pval,q2,mae_loo,var_ratio,noise_level,"
                "rho_matches_artifact\n")
        matches_artifact = int(abs(rho - control_rho) < 0.1)
        f.write(f"{args.state_label} {n} {var_y:.4f} {np.sqrt(var_y):.4f} "
                f"{y.max() - y.min():.4f} {control_rho:.4f} {control_mae:.4f} "
                f"{rho:.4f} {pval:.4f} {q2:.4f} {mae_loo:.4f} {var_ratio:.4f} "
                f"{noise_level:.4f} {matches_artifact}\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
