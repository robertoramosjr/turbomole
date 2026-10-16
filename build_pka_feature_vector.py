#!/usr/bin/env python3
"""
build_pka_feature_vector.py

Assembles the candidate descriptor vector for Artepilin C's pKa
calibration survey (protonation_and_conformer_investigation.md
Sec. 9.30) -- population-weighted average, over each state's LABELED
(real HSE06) subset, of every descriptor actually available in this
project right now, for each of the 3 acid/base pairs:
  - neutral <-> monoanion_a (carboxyl deprotonation)
  - neutral <-> monoanion_b (phenol deprotonation)
  - monoanion_a <-> dianion (second/phenol deprotonation, stepwise
    from the dominant carboxylate tautomer)

IMPORTANT -- this does NOT compute pKa. It only builds the descriptor
table and reports, explicitly, which of the originally-requested
descriptor types are and are NOT available without new calculations:

  - delta_G via `freeh`: NOT AVAILABLE -- `freeh` has never been run
    anywhere in this project (checked: no freeh output, no script
    calls it). Would need a new post-processing step on the existing
    aoforce Hessians (cheap relative to DFT, but still new work,
    outside this task's scope).
  - Bader charge / bond order on the IONIZABLE O-H: NOT AVAILABLE as
    such -- `population_feature.csv` (Sec. 9.10) only tracked the
    CARBONYL C/O (B1's target atom), which is a DIFFERENT atom than
    the ionizable carboxyl/phenol oxygen for most conformers. Reported
    below under its real name (carbonyl, not ionizable-site) rather
    than mislabeled.
  - Bond order O-H specifically: NOT AVAILABLE -- never computed for
    any O-H bond, ionizable or not.
  - HOMO-LUMO gap: AVAILABLE (extract_homo_lumo_gap.py, parsed for
    free from each labeled conformer's existing ridft.out -- no new
    calculation).
  - COSMO moments: NOT AVAILABLE -- COSMO/implicit solvation was never
    used anywhere in this project; every calculation so far is gas-phase.

Usage:
    python ~/work_turbomole/scripts/build_pka_feature_vector.py \
        --state-dir neutral=artepillin_c_neutral_d3_hse06 \
                    monoanion_a=artepillin_c_monoanion_a_d3_hse06 \
                    monoanion_b=artepillin_c_monoanion_b_d3_hse06 \
                    dianion=artepillin_c_dianion_d3_hse06 \
        --output pka_feature_vector.csv
"""

import argparse
import csv
import os


def load_population_weighted_avg(state_dir, labeled_dirs):
    pop_path = os.path.join(state_dir, "layer_a", "population_table_dft.csv")
    with open(pop_path) as f:
        pop = {r["job_dir"]: float(r["r2scan3c_population"]) for r in csv.DictReader(f)}

    gap_path = os.path.join(state_dir, "layer_b", "homo_lumo_gap.csv")
    with open(gap_path) as f:
        gap = {r["job_dir"]: float(r["homo_lumo_gap_eV"]) for r in csv.DictReader(f)}

    pop_feat_path = os.path.join(state_dir, "layer_b", "population_feature.csv")
    carbonyl_bader = {}
    carbonyl_bond_order = {}
    if os.path.isfile(pop_feat_path):
        with open(pop_feat_path) as f:
            for r in csv.DictReader(f):
                carbonyl_bader[r["job_dir"]] = float(r["bader_charge_carbonyl_o"])
                carbonyl_bond_order[r["job_dir"]] = float(r["co_bond_order_sen"])

    common = [jd for jd in labeled_dirs if jd in pop and jd in gap]
    w_sum = sum(pop[jd] for jd in common)

    def w_avg(d):
        vals = [(jd, d[jd]) for jd in common if jd in d]
        if not vals:
            return None
        w = sum(pop[jd] for jd, _ in vals)
        return sum(pop[jd] * v for jd, v in vals) / w if w > 0 else None

    return {
        "n_labeled": len(common),
        "population_covered": w_sum,
        "homo_lumo_gap_eV": w_avg(gap),
        "carbonyl_bader_charge": w_avg(carbonyl_bader),
        "carbonyl_bond_order_sen": w_avg(carbonyl_bond_order),
        "delta_g_freeh_eV": None,  # not available -- see docstring
        "ionizable_oh_bond_order": None,  # not available
        "cosmo_dipole": None,  # not available
    }


def main():
    parser = argparse.ArgumentParser(
        description="Assemble the pKa candidate descriptor vector from "
                    "already-available data, flagging what is missing.")
    parser.add_argument("--state-dir", nargs="+", required=True,
                         help="state=directory pairs, e.g. "
                              "neutral=artepillin_c_neutral_d3_hse06")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    state_dirs = {}
    for entry in args.state_dir:
        state, path = entry.split("=", 1)
        state_dirs[state] = path

    descriptors = {}
    for state, path in state_dirs.items():
        labeled_dirs = os.listdir(os.path.join(path, "layer_b", "labeled_jobs"))
        descriptors[state] = load_population_weighted_avg(path, labeled_dirs)
        print(f"=== {state} (N_labeled={descriptors[state]['n_labeled']}, "
              f"population covered={descriptors[state]['population_covered']:.4f}) ===")
        for k, v in descriptors[state].items():
            if k in ("n_labeled", "population_covered"):
                continue
            print(f"  {k}: {v if v is not None else 'NOT AVAILABLE'}")
        print()

    pairs = [
        ("neutral", "monoanion_a", "carboxyl deprotonation"),
        ("neutral", "monoanion_b", "phenol deprotonation"),
        ("monoanion_a", "dianion", "second (phenol) deprotonation, stepwise from carboxylate"),
    ]

    fields = ["homo_lumo_gap_eV", "carbonyl_bader_charge", "carbonyl_bond_order_sen",
              "delta_g_freeh_eV", "ionizable_oh_bond_order", "cosmo_dipole"]

    rows = []
    print("=== Pair deltas (base - acid) ===")
    for acid, base, label in pairs:
        if acid not in descriptors or base not in descriptors:
            continue
        row = {"pair": f"{acid}->{base}", "description": label}
        for field in fields:
            a, b = descriptors[acid][field], descriptors[base][field]
            row[field] = (b - a) if (a is not None and b is not None) else ""
        rows.append(row)
        print(f"  {acid} -> {base} ({label}):")
        for field in fields:
            print(f"    delta_{field}: {row[field] if row[field] != '' else 'NOT AVAILABLE'}")

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["pair", "description"] + fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {args.output}")
    print("\n*** pKa was NOT computed -- this is descriptor assembly only. ***")


if __name__ == "__main__":
    main()
