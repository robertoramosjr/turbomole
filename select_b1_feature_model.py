#!/usr/bin/env python3
"""
select_b1_feature_model.py

Automates the feature/model selection protocol agreed for B1 after the
monoanion B audit (Sec. 9.8-9.9 of protonation_and_conformer_investigation.md),
so it is applied consistently to the neutral/monoanion A feature
expansion instead of re-deriving the decision by hand each time:

  1. Pearson correlation of EACH candidate feature IN ISOLATION against
     the delta target. Report r and p-value for every feature, not
     just the ones that look promising -- a feature with no individual
     correlation is not a reason on its own to skip it from step 2 (a
     purely nonlinear relationship could still be missed by Pearson),
     but it lowers the prior that a linear model will use it well.
  2. Every feature gets an OLS (1 feature, 2 parameters) leave-one-out
     evaluation on the SAME fold structure as the rest of this
     project's LOO analyses (diagnose_b1_surrogate.py /
     train_b1_gpr.py / test_b1_linear_baseline.py). Q2 > 0 is the bar
     for "this feature alone has real learnable signal."
  3. If two or more features individually clear Q2 > 0, a combined
     linear (no interaction term) LOO-OLS is also tried. The combo is
     recommended ONLY if its Q2 beats the best single feature's Q2 --
     otherwise the best single feature wins and the combo is reported
     but explicitly marked "do not use" (the same rule that closed the
     monoanion B model: 2 features scored worse than 1 there).
  4. A GPR (ARD Matern 5/2 + WhiteKernel, matching
     diagnose_b1_surrogate.py) is recommended ONLY if EITHER:
       (a) no linear model (single or combined) reaches Q2 > 0, AND
           there is a concrete reason to suspect nonlinearity (this
           script does not decide that automatically -- it prints the
           gate and requires a human judgment call, logged in the
           output), OR
       (b) N is large enough that a GPR's extra hyperparameters
           (2 length scales + signal variance + noise level for a
           2-feature kernel) have a reasonable points-per-parameter
           ratio -- this script uses a configurable rule of thumb
           (--min-n-per-gpr-param, default 5) and reports whether N
           clears it, but does not silently run a GPR search itself.
     A multi-feature GPR is never the first thing this script tries.

Usage:
    python ~/work_turbomole/scripts/select_b1_feature_model.py \
        --training-data layer_b/b1_training_data.csv \
        --features layer_b/oh_distance_feature.csv \
                   layer_b/planarity_feature.csv \
                   layer_b/bond_order_feature.csv \
                   layer_b/bader_charge_feature.csv \
        --state-label neutral \
        --output layer_b/b1_feature_selection.dat
"""

import argparse
import csv

import numpy as np
from scipy.stats import pearsonr


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


def loo_ols(X, y):
    n = len(y)
    preds = np.zeros(n)
    for i in range(n):
        train_idx = [j for j in range(n) if j != i]
        design = np.column_stack([np.ones(len(train_idx)), X[train_idx]])
        coef, *_ = np.linalg.lstsq(design, y[train_idx], rcond=None)
        x_test = np.concatenate([[1.0], X[i]])
        preds[i] = x_test @ coef
    return preds


def q2(y, preds):
    ss_tot = np.sum((y - y.mean()) ** 2)
    ss_res = np.sum((y - preds) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def main():
    parser = argparse.ArgumentParser(
        description="Apply the B1 feature/model selection protocol "
                    "(Pearson -> single-feature LOO-OLS -> combo only if "
                    "it beats the best single -> GPR only if gated conditions "
                    "hold) to a new set of candidate features.")
    parser.add_argument("--training-data", required=True)
    parser.add_argument("--features", required=True, nargs="+")
    parser.add_argument("--target-column", default="delta_cm1",
                         help="Column name in --training-data holding the "
                              "delta target [default: delta_cm1, B1's "
                              "carbonyl shift -- pass e.g. delta_energy_eV "
                              "or delta_log10_f for B2]")
    parser.add_argument("--state-label", required=True)
    parser.add_argument("--min-n-per-gpr-param", type=int, default=5,
                         help="Rule-of-thumb minimum labeled points per GPR "
                              "hyperparameter before a multi-feature GPR is "
                              "considered reasonable [default: 5]")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    job_dirs, X, y, feature_names = load_data(
        args.training_data, args.features, args.target_column)
    n = len(job_dirs)
    print(f"=== {args.state_label} (N={n}) ===\n")

    if n < 4:
        print("N < 4: too few points for this protocol.")
        with open(args.output, "w") as f:
            f.write(f"# {args.state_label}: N={n}, skipped\n")
        return

    # Step 1 + 2: every feature, isolated
    print("--- Step 1+2: per-feature Pearson correlation and LOO-OLS ---")
    single_q2 = {}
    single_results = []
    for i, name in enumerate(feature_names):
        x = X[:, i]
        r, p = pearsonr(x, y)
        preds = loo_ols(x.reshape(-1, 1), y)
        q2_val = q2(y, preds)
        single_q2[name] = q2_val
        flag = "Q2>0 -- candidate" if q2_val > 0 else "Q2<=0 -- no linear signal alone"
        print(f"  {name:30s} r={r:+.4f} (p={p:.4f})  Q2_loo_ols={q2_val:+.4f}  [{flag}]")
        single_results.append((name, r, p, q2_val))

    candidates = [name for name, q in single_q2.items() if q > 0]
    best_single = max(single_q2, key=single_q2.get) if single_q2 else None
    best_single_q2 = single_q2[best_single] if best_single else float("-inf")

    # Step 3: combo, only reported (and only recommended) if it beats best single
    combo_q2 = None
    if len(candidates) >= 2:
        print(f"\n--- Step 3: {len(candidates)} features cleared Q2>0 "
              f"({candidates}) -- trying the linear combination ---")
        idx = [feature_names.index(c) for c in candidates]
        preds_combo = loo_ols(X[:, idx], y)
        combo_q2 = q2(y, preds_combo)
        print(f"  Combined ({len(candidates)} features, linear, no interaction): "
              f"Q2 = {combo_q2:+.4f}")
        if combo_q2 > best_single_q2:
            print(f"  => Combo BEATS best single feature ({best_single}, "
                  f"Q2={best_single_q2:+.4f}) -- recommend the combo.")
        else:
            print(f"  => Combo does NOT beat best single feature ({best_single}, "
                  f"Q2={best_single_q2:+.4f}) -- recommend {best_single} alone, "
                  f"discard the combo (same rule that closed monoanion B).")
    elif len(candidates) == 1:
        print(f"\n--- Step 3: only one feature ({candidates[0]}) cleared Q2>0 "
              f"-- nothing to combine, it is the recommendation. ---")
    else:
        print("\n--- Step 3: no feature cleared Q2>0 in isolation -- "
              "no combination attempted (combining non-signal features "
              "does not create signal). ---")

    # Step 4: GPR gating -- reported, never auto-run
    n_features_for_gpr = len(feature_names)
    n_gpr_params = 2 * n_features_for_gpr + 2  # length scales + signal var + noise
    min_n_needed = args.min_n_per_gpr_param * n_gpr_params
    print(f"\n--- Step 4: GPR gating (informational only -- this script does "
          f"not run a GPR) ---")
    print(f"  A {n_features_for_gpr}-feature ARD GPR has ~{n_gpr_params} "
          f"hyperparameters; rule of thumb wants N >= "
          f"{args.min_n_per_gpr_param} x that = {min_n_needed}.")
    print(f"  Current N = {n} -- "
          f"{'clears' if n >= min_n_needed else 'does NOT clear'} that bar.")
    if candidates:
        print(f"  A linear model already has Q2>0 ({best_single} or the combo) "
              f"-- per this protocol, do not reach for a GPR unless there is a "
              f"specific, separately-documented reason to suspect nonlinearity "
              f"the linear fit is missing.")
    else:
        print(f"  No linear model reached Q2>0 -- a GPR is only justified here "
              f"if (a) there is a concrete reason to suspect nonlinearity "
              f"(human judgment call, not automatic), and ideally (b) N clears "
              f"the rule-of-thumb bar above. Do not default to a multi-feature "
              f"GPR just because linear models failed -- that repeats the "
              f"monoanion B mistake of over-parameterizing a small sample.")

    with open(args.output, "w") as f:
        f.write("descriptor feature,pearson_r,pearson_p,q2_loo_ols,is_candidate\n")
        for name, r, p, qv in single_results:
            f.write(f"{name} {r:.4f} {p:.4f} {qv:.4f} {int(qv > 0)}\n")
        if combo_q2 is not None:
            f.write(f"combined_linear NA NA {combo_q2:.4f} "
                     f"{int(combo_q2 > best_single_q2)}\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
