#!/usr/bin/env python3
"""
select_random_population_weighted_batch.py

Draws a confirmation batch by population-weighted RANDOM sampling
(without replacement) from a pool of unconfirmed conformers, instead
of by extreme predicted value -- the control needed to separate real
offset error from winner's-curse/extreme-order-statistics selection
bias (protonation_and_conformer_investigation.md Sec. 9.24, Task 4).
Sampling weight is each conformer's own r2SCAN-3c Boltzmann population
(keeps the sample relevant to the final spectrum -- a random UNWEIGHTED
draw would waste slots on near-zero-population conformers), but
NEVER the conformer's own predicted energy/f, which is exactly the
quantity being validated.

Usage:
    python ~/work_turbomole/scripts/select_random_population_weighted_batch.py \
        --contribution-rank layer_b/b2_spectral_contribution_rank.csv \
        --exclude-already-confirmed job_0018 job_0010 job_0023 ... \
        --pool-method surrogate_max_f --pool-uncertain-inflated False \
        --n-select 8 --seed 0 \
        --output layer_b/b2_random_control_batch.csv
"""

import argparse
import csv

import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="Population-weighted random (not extreme-value) "
                    "confirmation batch draw, to isolate winner's-curse "
                    "selection bias from real offset error.")
    parser.add_argument("--contribution-rank", required=True,
                         help="rank_b2_spectral_contribution.py output")
    parser.add_argument("--exclude-already-confirmed", nargs="*", default=[])
    parser.add_argument("--pool-method", default="surrogate_max_f",
                         help="Only draw from rows with this 'method' value")
    parser.add_argument("--pool-uncertain-inflated", default="False",
                         help="Only draw from rows with this 'uncertain_inflated' "
                              "value (as printed in the CSV, e.g. 'False')")
    parser.add_argument("--n-select", type=int, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    excluded = set(args.exclude_already_confirmed)
    with open(args.contribution_rank) as f:
        rows = [r for r in csv.DictReader(f)
                if r["job_dir"] not in excluded
                and r["method"] == args.pool_method
                and r["uncertain_inflated"] == args.pool_uncertain_inflated]

    if len(rows) < args.n_select:
        raise SystemExit(f"Pool has only {len(rows)} eligible conformers, "
                          f"cannot select {args.n_select}.")

    weights = np.array([float(r["population"]) for r in rows])
    weights = weights / weights.sum()

    rng = np.random.default_rng(args.seed)
    idx = rng.choice(len(rows), size=args.n_select, replace=False, p=weights)
    selected = [rows[i] for i in idx]

    print(f"Pool size: {len(rows)} (method={args.pool_method}, "
          f"uncertain_inflated={args.pool_uncertain_inflated})")
    print(f"Drew {args.n_select} conformers, population-weighted, "
          f"WITHOUT using their own predicted value:")
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["job_dir", "population", "predicted_energy_eV",
                           "predicted_osc", "contribution"])
        writer.writeheader()
        for r in selected:
            writer.writerow({
                "job_dir": r["job_dir"], "population": r["population"],
                "predicted_energy_eV": r["predicted_energy_eV"],
                "predicted_osc": r["predicted_osc"],
                "contribution": r["contribution"],
            })
            print(f"  {r['job_dir']}: pop={float(r['population']):.4f}, "
                  f"predicted_f={float(r['predicted_osc']):.4f}")

    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
