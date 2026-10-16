#!/usr/bin/env python3
"""
simulate_offset_selection_bias.py

Checks how much of an observed positive bias in a small confirmation
batch could be pure "winner's curse" / extreme-order-statistics
artifact from selecting conformers by their PREDICTED value (energy or
f, which both ride on the same flat, supposedly zero-mean-residual
offset -- protonation_and_conformer_investigation.md Sec. 9.21/9.23),
rather than evidence of real, additional bias in the offset itself.

Model: --pool-size conformers each have an unknown true residual
(HSE06 minus offset-corrected-sTDA) drawn i.i.d. from
N(0, --residual-std) -- the offset's OWN LOO-calibrated std, i.e. the
null hypothesis is "the offset is honestly unbiased and homoscedastic,
and the batch just happened to draw an extreme top-K by chance because
selection targeted extreme predicted values". This is the WORST-CASE
version of that hypothesis: it assumes selection order tracks the
residual perfectly (selecting the top predicted magnitude means
selecting the top residual too), which is more favorable to the
"it's just selection" explanation than reality, where selection was by
predicted f/contribution, only loosely coupled to the energy residual.
If even this worst-case selection model cannot reproduce the observed
bias, that is strong evidence of a real, additional bias source beyond
pure order-statistics selection.

For each of --n-trials repeats: draw --pool-size samples, take the
mean of the top --n-selected (largest) values, record it. Reports
where --observed-mean-bias falls in that simulated null distribution.

Usage:
    python ~/work_turbomole/scripts/simulate_offset_selection_bias.py \
        --pool-size 118 --n-selected 5 --residual-std 0.2140 \
        --observed-mean-bias 0.0985 \
        --n-trials 5000 --seed 0 \
        --output layer_b/selection_bias_simulation_population_weighted.dat
"""

import argparse

import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="Simulate the worst-case winner's-curse bias from "
                    "selecting the top-K conformers by predicted value out "
                    "of a pool with i.i.d. zero-mean Gaussian residuals, "
                    "and compare against an observed mean bias.")
    parser.add_argument("--pool-size", type=int, required=True)
    parser.add_argument("--n-selected", type=int, required=True)
    parser.add_argument("--residual-std", type=float, required=True,
                         help="The offset model's own LOO/calibration std "
                              "(the assumed per-conformer residual std "
                              "under the null of 'honestly unbiased')")
    parser.add_argument("--observed-mean-bias", type=float, required=True)
    parser.add_argument("--n-trials", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.n_trials < 1000:
        raise SystemExit("--n-trials must be at least 1000 for a stable tail estimate.")

    rng = np.random.default_rng(args.seed)
    top_k_means = np.empty(args.n_trials)
    for t in range(args.n_trials):
        draws = rng.normal(0.0, args.residual_std, size=args.pool_size)
        top_k = np.sort(draws)[-args.n_selected:]
        top_k_means[t] = top_k.mean()

    percentile = float(np.mean(top_k_means <= args.observed_mean_bias) * 100)
    p_value = float(np.mean(top_k_means >= args.observed_mean_bias))

    print(f"Pool size={args.pool_size}, n_selected={args.n_selected}, "
          f"residual_std={args.residual_std}, n_trials={args.n_trials}")
    print(f"Worst-case selection-only null: mean(top-{args.n_selected} of "
          f"{args.pool_size}) = {top_k_means.mean():.4f} +/- "
          f"{top_k_means.std():.4f} (std across trials)")
    print(f"Null distribution: p5={np.percentile(top_k_means, 5):.4f}, "
          f"median={np.median(top_k_means):.4f}, "
          f"p95={np.percentile(top_k_means, 95):.4f}, "
          f"max={top_k_means.max():.4f}")
    print(f"\nObserved mean bias = {args.observed_mean_bias:.4f}")
    print(f"  sits at percentile {percentile:.1f} of the worst-case "
          f"selection-only null (fraction of null >= observed = {p_value:.4f})")

    if args.observed_mean_bias > np.percentile(top_k_means, 95):
        verdict = ("Observed bias EXCEEDS the 95th percentile of even the "
                   "worst-case (selection order = residual order) null -- "
                   "pure selection-effect cannot plausibly explain it; "
                   "this points to a real, additional bias source.")
    elif args.observed_mean_bias < 0:
        verdict = "Observed bias is not even positive -- no selection-effect question to resolve."
    else:
        verdict = ("Observed bias falls within (or below) the range the "
                   "worst-case selection-only null can produce -- "
                   "selection effect alone is a plausible explanation, "
                   "though this does not prove it is the ONLY explanation "
                   "(real bias and selection effect are not mutually exclusive).")
    print(f"\n=> {verdict}")

    with open(args.output, "w") as f:
        f.write("descriptor pool_size,n_selected,residual_std,n_trials,"
                "null_mean,null_std,null_p5,null_median,null_p95,"
                "observed_mean_bias,observed_percentile,empirical_pvalue\n")
        f.write(f"{args.pool_size} {args.n_selected} {args.residual_std} "
                f"{args.n_trials} {top_k_means.mean():.4f} {top_k_means.std():.4f} "
                f"{np.percentile(top_k_means, 5):.4f} {np.median(top_k_means):.4f} "
                f"{np.percentile(top_k_means, 95):.4f} {args.observed_mean_bias:.4f} "
                f"{percentile:.2f} {p_value:.4f}\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
