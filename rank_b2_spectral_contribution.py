#!/usr/bin/env python3
"""
rank_b2_spectral_contribution.py

Ranks conformers by how much each one contributes to the integrated
S1 band intensity of a protonation state's population-weighted UV-Vis
spectrum -- population_weight * oscillator_strength -- using the
already-established r2SCAN-3c Boltzmann population
(layer_a/population_table_dft.csv, the same "r2scan3c_population"
column build_population_weighted_spectrum.py already uses for
production spectra) and B2's per-conformer S1 prediction
(apply_b2_osc_surrogate.py's output: real HSE06 label, forced-state
override, or offset-corrected surrogate, whichever applies).

This is NOT the final spectrum -- it is a ranking of "how much each
conformer's own confirmation/error would move the aggregate result",
used to prioritize which unconfirmed conformers are worth a real
HSE06/TD-DFT calculation before the spectrum is finalized.
Conformers marked 'excluded' in the B2 predictions table (low sTDA
selection confidence, no manual resolution yet) have no usable f and
are reported separately, not silently dropped from view.

Usage:
    python ~/work_turbomole/scripts/rank_b2_spectral_contribution.py \
        --population-table layer_a/population_table_dft.csv \
        --b2-predictions layer_b/b2_osc_predictions_full_ensemble.csv \
        --state-label neutral \
        --top-n 20 \
        --output layer_b/b2_spectral_contribution_rank.csv
"""

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(
        description="Rank conformers by population * oscillator-strength "
                    "contribution to a protonation state's aggregate S1 band.")
    parser.add_argument("--population-table", required=True)
    parser.add_argument("--b2-predictions", required=True)
    parser.add_argument("--state-label", required=True)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.population_table) as f:
        population = {row["job_dir"]: float(row["r2scan3c_population"])
                      for row in csv.DictReader(f)}

    with open(args.b2_predictions) as f:
        predictions = {row["job_dir"]: row for row in csv.DictReader(f)}

    rows = []
    excluded = []
    for jd, pop in population.items():
        pred = predictions.get(jd)
        if pred is None:
            print(f"WARNING: {jd} has a population weight but no B2 "
                  f"prediction row, skipping")
            continue
        if pred["excluded"] == "True" or not pred["predicted_osc"]:
            excluded.append((jd, pop))
            continue
        f_val = float(pred["predicted_osc"])
        e_val = float(pred["predicted_energy_eV"])
        contribution = pop * f_val
        rows.append({
            "job_dir": jd, "population": pop, "predicted_energy_eV": e_val,
            "predicted_osc": f_val, "contribution": contribution,
            "method": pred["method"], "uncertain_inflated": pred["uncertain_inflated"],
        })

    rows.sort(key=lambda r: r["contribution"], reverse=True)
    total_contribution = sum(r["contribution"] for r in rows)
    total_excluded_population = sum(p for _, p in excluded)

    print(f"=== {args.state_label} spectral contribution ranking (N={len(rows)}"
          f" ranked, {len(excluded)} excluded from ranking) ===")
    print(f"Total excluded population (no usable f): "
          f"{total_excluded_population:.4f} "
          f"({100 * total_excluded_population:.1f}% of the state's total weight)")
    print(f"\nTop-{args.top_n} by contribution (population x f):")
    for r in rows[:args.top_n]:
        pct = 100 * r["contribution"] / total_contribution if total_contribution else 0
        print(f"  {r['job_dir']}: pop={r['population']:.4f}, "
              f"f={r['predicted_osc']:.4f}, E={r['predicted_energy_eV']:.3f} eV, "
              f"contribution={r['contribution']:.5f} ({pct:.1f}% of total), "
              f"method={r['method']}"
              + (" [UNCERTAIN_INFLATED]" if r["uncertain_inflated"] == "True" else ""))

    with open(args.output, "w", newline="") as f:
        fieldnames = ["job_dir", "population", "predicted_energy_eV",
                      "predicted_osc", "contribution", "method", "uncertain_inflated"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {args.output} ({len(rows)} ranked rows)")
    if excluded:
        excluded.sort(key=lambda x: x[1], reverse=True)
        print(f"Excluded conformers (population, sorted): "
              + ", ".join(f"{jd}({p:.4f})" for jd, p in excluded[:10])
              + (" ..." if len(excluded) > 10 else ""))


if __name__ == "__main__":
    main()
