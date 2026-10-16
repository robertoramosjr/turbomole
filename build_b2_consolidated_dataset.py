#!/usr/bin/env python3
"""
build_b2_consolidated_dataset.py

Builds one consolidated table, one row per conformer, across all 4
protonation states, with everything needed to build the final B2
UV-Vis figures: S1 energy + oscillator strength (real HSE06 where
available, corrected surrogate/offset otherwise -- exactly as decided
and documented per state in protonation_and_conformer_investigation.md
Sec. 9.15-9.28), and the r2SCAN-3c Boltzmann population weight
(population_table_dft.csv, already used throughout Layer A/B).

Per-state source of (E, f):
  - neutral: apply_b2_osc_surrogate.py v2 output (real HSE06, forced-
    state, or slope-intercept surrogate, per conformer).
  - monoanion_b: real HSE06 for all 11 conformers (100% of that
    state's ensemble is already labeled -- no surrogate needed).
  - monoanion_a: real HSE06 for the 20 labeled conformers only (covers
    26.6%% of the state's population) -- NO validated B2 model exists
    for the other 123 (Sec. 9.18: no feature survived calibration for
    either target in this state), so those rows are left OUT rather
    than filled with an unvalidated guess. Flagged explicitly as
    partial coverage.
  - dianion: real HSE06 for all 3 conformers (100% of that state's
    ensemble -- N=3 is too small for any model, but happens to have
    complete real coverage anyway).

A "delta_g_eV" column is included in the header for forward
compatibility with the pH-speciation figure, but is left EMPTY for
every row: no pKa/free-energy-of-deprotonation calibration exists yet
anywhere in this project (checked -- no file matches "pka"/"delta_g"/
"speciation" under this state's directories, and the referenced
`plano_especiacao_ph_alizarina.md` planning document is not present on
this machine). Populating it is a separate, not-yet-started prerequisite
(see MEMORY.md's pending-tasks note on pKa calibration), not something
this script can compute from what already exists.

Usage:
    python ~/work_turbomole/scripts/build_b2_consolidated_dataset.py \
        --neutral-predictions artepillin_c_neutral_d3_hse06/layer_b/b2_osc_predictions_full_ensemble_v2.csv \
        --neutral-population artepillin_c_neutral_d3_hse06/layer_a/population_table_dft.csv \
        --monoanion-b-hse06 artepillin_c_monoanion_b_d3_hse06/layer_b/s1_hse06.csv \
        --monoanion-b-population artepillin_c_monoanion_b_d3_hse06/layer_a/population_table_dft.csv \
        --monoanion-a-hse06 artepillin_c_monoanion_a_d3_hse06/layer_b/s1_hse06.csv \
        --monoanion-a-population artepillin_c_monoanion_a_d3_hse06/layer_a/population_table_dft.csv \
        --dianion-hse06 artepillin_c_dianion_d3_hse06/layer_b/s1_hse06.csv \
        --dianion-population artepillin_c_dianion_d3_hse06/layer_a/population_table_dft.csv \
        --output b2_consolidated_dataset.csv
"""

import argparse
import csv


def load_population(path):
    with open(path) as f:
        return {r["job_dir"]: float(r["r2scan3c_population"]) for r in csv.DictReader(f)}


def rows_from_v2_predictions(pred_path, population_path, state):
    pop = load_population(population_path)
    rows = []
    with open(pred_path) as f:
        for r in csv.DictReader(f):
            if r["excluded"] == "True" or not r["predicted_osc"]:
                continue
            jd = r["job_dir"]
            rows.append({
                "job_dir": jd, "protonation_state": state,
                "energy_eV": r["predicted_energy_eV"], "osc_f": r["predicted_osc"],
                "source": r["method"], "population": pop.get(jd, ""),
                "state_selection_mismatch_rate": r.get("state_selection_mismatch_rate", ""),
                "delta_g_eV": "",
            })
    return rows


def rows_from_real_hse06_only(hse06_path, population_path, state):
    pop = load_population(population_path)
    rows = []
    with open(hse06_path) as f:
        for r in csv.DictReader(f):
            jd = r["job_dir"]
            rows.append({
                "job_dir": jd, "protonation_state": state,
                "energy_eV": r["s1_energy_eV"], "osc_f": r["s1_osc_len"],
                "source": "hse06_ground_truth", "population": pop.get(jd, ""),
                "state_selection_mismatch_rate": "",
                "delta_g_eV": "",
            })
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Consolidate per-conformer S1 (E, f) + population across "
                    "all 4 protonation states into one table for the final "
                    "B2 UV-Vis figures.")
    parser.add_argument("--neutral-predictions", required=True)
    parser.add_argument("--neutral-population", required=True)
    parser.add_argument("--monoanion-b-hse06", required=True)
    parser.add_argument("--monoanion-b-population", required=True)
    parser.add_argument("--monoanion-a-hse06", required=True)
    parser.add_argument("--monoanion-a-population", required=True)
    parser.add_argument("--dianion-hse06", required=True)
    parser.add_argument("--dianion-population", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    all_rows = []
    all_rows += rows_from_v2_predictions(
        args.neutral_predictions, args.neutral_population, "neutral")
    all_rows += rows_from_real_hse06_only(
        args.monoanion_b_hse06, args.monoanion_b_population, "monoanion_b")
    all_rows += rows_from_real_hse06_only(
        args.monoanion_a_hse06, args.monoanion_a_population, "monoanion_a")
    all_rows += rows_from_real_hse06_only(
        args.dianion_hse06, args.dianion_population, "dianion")

    with open(args.output, "w", newline="") as f:
        fieldnames = ["job_dir", "protonation_state", "energy_eV", "osc_f",
                      "source", "population", "state_selection_mismatch_rate",
                      "delta_g_eV"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"Wrote {len(all_rows)} rows -> {args.output}")
    for state in ("neutral", "monoanion_a", "monoanion_b", "dianion"):
        n = sum(1 for r in all_rows if r["protonation_state"] == state)
        pop_sum = sum(float(r["population"]) for r in all_rows
                       if r["protonation_state"] == state and r["population"] != "")
        print(f"  {state}: {n} conformers, population coverage sum = {pop_sum:.4f}")
    print(f"\nNOTE: delta_g_eV is EMPTY for every row -- no pKa/speciation "
          f"data exists yet in this project (see script docstring). The "
          f"pH-speciation figure cannot be built until that is computed "
          f"separately.")
    print(f"NOTE: monoanion_a covers only its 20 labeled conformers "
          f"(26.6% of that state's population) -- no B2 model was "
          f"validated for the rest (Sec. 9.18), so they are left out "
          f"rather than filled with an unvalidated guess.")


if __name__ == "__main__":
    main()
