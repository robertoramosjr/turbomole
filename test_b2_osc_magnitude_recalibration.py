#!/usr/bin/env python3
"""
test_b2_osc_magnitude_recalibration.py

Tests whether B2's flat additive offset in log10(f) space
(delta_log10_f = log10(f_HSE06) - log10(f_sTDA), applied as a single
constant mean+std -- Sec. 9.21) is the right model, or whether the
true relationship has a slope != 1 (protonation_and_conformer_investigation.md
Sec. 9.25): fits log10(f_HSE06) = a + b*log10(f_sTDA) by OLS on a
CLEAN pool (character-matched sTDA on both sides, no state-selection
mismatches, no manual special cases, no duplicates), tests b=1 via the
coefficient's own confidence interval, and validates whichever model
wins (b=1, i.e. keep the flat offset, vs. b!=1, i.e. adopt the fitted
slope) with this project's standard LOO + permutation protocol
(feature_selection_protocol.md) -- never skipping the validation step
just because the fit came from a regression instead of a single
feature search.

Also applies the winning model out-of-sample to the extreme-value
confirmation batch (excluding job_0023, a character MISMATCH -- a
different error category, not a magnitude question -- and job_0018, a
duplicate of an already-labeled conformer) as an independent check:
does the new model reduce that batch's observed bias too, even though
those points were selected by extreme predicted value (enriched by
winner's-curse, Sec. 9.24 Task 1) rather than drawn cleanly?

Usage:
    python ~/work_turbomole/scripts/test_b2_osc_magnitude_recalibration.py \
        --labeled-delta-target layer_b/b2_osc_delta_target.csv \
        --random-control-delta-target layer_b/b2_osc_delta_target_random_control.csv \
        --extreme-batch-hse06 layer_b/s1_hse06_confirmation_batch.csv \
        --extreme-batch-stda-confidence layer_b/stda_confidence_full_ensemble.csv \
        --extreme-batch-exclude job_0023 job_0018 \
        --n-permutations 1000 --seed 0 \
        --output layer_b/b2_osc_magnitude_recalibration.dat
"""

import argparse
import csv
import math

import numpy as np


def loo_ols_two_param(x, y):
    """Exact LOO for a 1-feature OLS (intercept + slope), via the
    PRESS/hat-matrix shortcut -- same closed form used throughout this
    project (permutation_test_b1_features.py)."""
    n = len(y)
    design = np.column_stack([np.ones(n), x])
    hat = design @ np.linalg.pinv(design.T @ design) @ design.T
    leverage = np.diag(hat)
    resid = y - hat @ y
    denom = np.clip(1 - leverage, 1e-10, None)
    loo_resid = resid / denom
    ss_tot = np.sum((y - y.mean()) ** 2)
    ss_res = np.sum(loo_resid ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def ols_fit_with_se(x, y):
    n = len(y)
    design = np.column_stack([np.ones(n), x])
    coef, residuals, rank, sv = np.linalg.lstsq(design, y, rcond=None)
    y_pred = design @ coef
    resid = y - y_pred
    dof = n - 2
    sigma2 = np.sum(resid ** 2) / dof
    cov = sigma2 * np.linalg.inv(design.T @ design)
    se = np.sqrt(np.diag(cov))
    return coef, se, dof


def main():
    parser = argparse.ArgumentParser(
        description="Test log10(f) slope vs. flat offset for B2's "
                    "oscillator-strength surrogate, validate with LOO + "
                    "permutation, and check out-of-sample on the extreme batch.")
    parser.add_argument("--labeled-delta-target", required=True,
                         help="build_b2_osc_delta_target.py output for the "
                              "original 20 labeled neutral conformers")
    parser.add_argument("--random-control-delta-target", required=True,
                         help="build_b2_osc_delta_target.py output for the "
                              "8 random-control conformers (character-matched)")
    parser.add_argument("--extreme-batch-hse06", required=True)
    parser.add_argument("--extreme-batch-stda-character-matched", required=True,
                         help="batch_stda_s1_character_matched.py output for "
                              "the extreme batch, if available; falls back to "
                              "--extreme-batch-stda-confidence (max-f, not "
                              "character-matched) if not given")
    parser.add_argument("--extreme-batch-exclude", nargs="*", default=[])
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.n_permutations < 500:
        raise SystemExit("--n-permutations must be at least 500.")

    # --- Task 1: assemble the clean 28-point pool ---
    pool = []
    for path in (args.labeled_delta_target, args.random_control_delta_target):
        with open(path) as f:
            for row in csv.DictReader(f):
                pool.append((row["job_dir"], float(row["hse06_osc"]), float(row["stda_osc"])))

    print(f"=== Task 1: clean pool ===")
    print(f"Assembled {len(pool)} clean (f_HSE06, f_sTDA) pairs "
          f"(20 original labeled + {len(pool) - 20} random control)")

    log_hse = np.array([math.log10(max(f_h, 1e-4)) for _, f_h, _ in pool])
    log_stda = np.array([math.log10(max(f_s, 1e-4)) for _, _, f_s in pool])

    # --- Task 2: slope vs. flat offset ---
    print(f"\n=== Task 2: log10(f_HSE06) = a + b*log10(f_sTDA) ===")
    coef, se, dof = ols_fit_with_se(log_stda, log_hse)
    a, b = coef
    se_a, se_b = se
    t_stat_b1 = (b - 1.0) / se_b
    from scipy import stats
    p_value_b1 = 2 * (1 - stats.t.cdf(abs(t_stat_b1), dof))
    ci95_b = (b - 1.96 * se_b, b + 1.96 * se_b)

    print(f"  a (intercept) = {a:+.4f} +/- {se_a:.4f}")
    print(f"  b (slope)     = {b:+.4f} +/- {se_b:.4f}  (95% CI [{ci95_b[0]:.4f}, {ci95_b[1]:.4f}])")
    print(f"  Test b=1: t={t_stat_b1:+.3f}, dof={dof}, p={p_value_b1:.4f}")

    slope_differs = p_value_b1 < 0.05
    if slope_differs:
        print(f"  => b differs significantly from 1 (p<0.05) -- adopting the "
              f"fitted slope+intercept model.")
        model_x, model_name = log_stda, "slope_and_intercept"
    else:
        print(f"  => b is NOT significantly different from 1 (p>=0.05) -- the "
              f"problem is the INTERCEPT, not the slope. Refitting a flat "
              f"offset (slope fixed at 1) on the clean 28-point pool.")
        model_x, model_name = None, "flat_offset_refit"

    if model_name == "flat_offset_refit":
        flat_offset_mean = float(np.mean(log_hse - log_stda))
        flat_offset_std = float(np.std(log_hse - log_stda, ddof=1))
        print(f"  Refit flat offset (N=28): mean={flat_offset_mean:+.4f}, "
              f"std={flat_offset_std:.4f} (compare to the original N=20 "
              f"offset used in Sec. 9.21)")
        y_pred_full = log_stda + flat_offset_mean

        def predict(x_stda_log):
            return x_stda_log + flat_offset_mean
    else:
        y_pred_full = a + b * log_stda

        def predict(x_stda_log):
            return a + b * x_stda_log

    # --- Task 3: LOO + permutation validation ---
    print(f"\n=== Task 3: validation (LOO + permutation, standard protocol) ===")
    if model_name == "slope_and_intercept":
        q2_loo = loo_ols_two_param(log_stda, log_hse)
    else:
        # LOO for a flat-offset (slope fixed =1) model: refit the mean
        # offset leaving each point out, exact closed form (mean of the
        # others), matches the LOO logic used for the neutral offset itself.
        resid_all = log_hse - log_stda
        n = len(resid_all)
        loo_pred = np.array([
            (resid_all.sum() - resid_all[i]) / (n - 1) + log_stda[i]
            for i in range(n)
        ])
        ss_tot = np.sum((log_hse - log_hse.mean()) ** 2)
        ss_res = np.sum((log_hse - loo_pred) ** 2)
        q2_loo = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    print(f"  Q2_LOO ({model_name}) = {q2_loo:+.4f}")

    rng = np.random.default_rng(args.seed)
    n = len(log_hse)
    null_q2 = np.zeros(args.n_permutations)
    for t in range(args.n_permutations):
        y_shuffled = rng.permutation(log_hse)
        if model_name == "slope_and_intercept":
            null_q2[t] = loo_ols_two_param(log_stda, y_shuffled)
        else:
            resid_s = y_shuffled - log_stda
            loo_pred_s = np.array([
                (resid_s.sum() - resid_s[i]) / (n - 1) + log_stda[i]
                for i in range(n)
            ])
            ss_tot_s = np.sum((y_shuffled - y_shuffled.mean()) ** 2)
            ss_res_s = np.sum((y_shuffled - loo_pred_s) ** 2)
            null_q2[t] = 1 - ss_res_s / ss_tot_s if ss_tot_s > 0 else float("nan")

    percentile = float(np.mean(null_q2 <= q2_loo) * 100)
    p_value = float(np.mean(null_q2 >= q2_loo))
    print(f"  Permutation ({args.n_permutations} shuffles): observed Q2 sits at "
          f"percentile {percentile:.1f}, empirical p={p_value:.4f}")
    if q2_loo > 0 and percentile >= 95:
        print(f"  => CONFIRMED: signal not explained by chance.")
    elif percentile < 10:
        print(f"  => WARNING: worse than nearly all permutations, investigate.")
    else:
        print(f"  => Consistent with the null -- treat with caution "
              f"(percentile {percentile:.1f} is not a confirmed signal).")

    # --- Out-of-sample check on the extreme batch ---
    print(f"\n=== Out-of-sample check: extreme-value confirmation batch ===")
    with open(args.extreme_batch_hse06) as f:
        extreme_hse = {r["job_dir"]: float(r["s1_osc_len"]) for r in csv.DictReader(f)}
    with open(args.extreme_batch_stda_character_matched) as f:
        extreme_stda = {r["job_dir"]: float(r["s1_osc_stda"]) for r in csv.DictReader(f)}

    excluded = set(args.extreme_batch_exclude)
    common = sorted((set(extreme_hse) & set(extreme_stda)) - excluded)
    print(f"  N={len(common)} (excluded: {sorted(excluded)})")

    old_offset_mean = 0.2066  # Sec. 9.21's original N=20 offset
    print(f"  {'job_dir':10s} {'f_real':>7s} {'old_pred':>9s} {'old_err':>8s} "
          f"{'new_pred':>9s} {'new_err':>8s}")
    old_abs_errs, new_abs_errs = [], []
    for jd in common:
        f_h = extreme_hse[jd]
        f_s = extreme_stda[jd]
        log_s = math.log10(max(f_s, 1e-4))
        old_pred_f = 10 ** (log_s + old_offset_mean)
        new_pred_f = 10 ** predict(log_s)
        old_err = old_pred_f - f_h
        new_err = new_pred_f - f_h
        old_abs_errs.append(abs(old_err))
        new_abs_errs.append(abs(new_err))
        print(f"  {jd:10s} {f_h:7.4f} {old_pred_f:9.4f} {old_err:+8.4f} "
              f"{new_pred_f:9.4f} {new_err:+8.4f}")

    print(f"\n  Mean abs error: old model = {np.mean(old_abs_errs):.4f}, "
          f"new model = {np.mean(new_abs_errs):.4f} "
          f"({'IMPROVED' if np.mean(new_abs_errs) < np.mean(old_abs_errs) else 'DID NOT improve'})")

    with open(args.output, "w") as f:
        f.write(f"descriptor model,n_pool,a,se_a,b,se_b,p_value_b_eq_1,"
                f"q2_loo,perm_percentile,perm_pvalue,"
                f"old_model_mean_abs_err_oos,new_model_mean_abs_err_oos\n")
        f.write(f"{model_name} {len(pool)} {a:.4f} {se_a:.4f} {b:.4f} {se_b:.4f} "
                f"{p_value_b1:.4f} {q2_loo:.4f} {percentile:.2f} {p_value:.4f} "
                f"{np.mean(old_abs_errs):.4f} {np.mean(new_abs_errs):.4f}\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
