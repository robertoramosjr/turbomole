#!/usr/bin/env python3
"""
rerank_layer_a_dft.py

Second half of Layer A: pulls the final r2SCAN-3c/def2-mTZVP electronic
energy out of each conformer job directory (produced by
fanout_conformer_jobs.py + job_r2scan3c_optimize.sh), and recomputes
Boltzmann populations *within the already GFN2-xTB-population-filtered
set* using these DFT energies instead. This is a re-ranking of Layer
A's survivors, not a new population filter -- conformers already
dropped by filter_boltzmann_population.py's 95% GFN2-xTB cutoff are
not reconsidered here.

Each conformer's Boltzmann weight is still multiplied by its rotamer
degeneracy (the "degeneracy" column filter_boltzmann_population.py
carried over from CREST's cre_members) -- same reasoning as there:
these are CREGEN-clustered near-identical rotamers of one conformer
family, not independent structures, so dropping the degeneracy weight
here would silently reintroduce the same under-weighting bug that was
fixed for the GFN2-xTB filter, just at the DFT re-ranking step instead.

Caveat (documented, not silently glossed over): this uses bare
electronic energy, not a free energy -- no vibrational/thermal
correction, since no Hessian has been computed at this cheap level.
That is an approximation appropriate for a screening-stage ranking,
not for the final production population (which should eventually use
a proper free energy once frequencies are available).

Usage:
    python ~/work_turbomole/scripts/rerank_layer_a_dft.py \
        --population-table layer_a/population_table.csv \
        --conformer-jobs layer_a/conformer_jobs \
        --output layer_a/population_table_dft.csv [--temp 298.15]
"""

import argparse
import csv
import os

KB_KCAL_PER_MOL_K = 0.0019872041
HARTREE_TO_KCAL = 627.5094740631


def read_last_energy(job_dir):
    energy_path = os.path.join(job_dir, "energy")
    if not os.path.isfile(energy_path):
        return None
    last_value = None
    with open(energy_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                last_value = float(parts[1])
    return last_value


def main():
    parser = argparse.ArgumentParser(
        description="Re-rank Layer A survivors by r2SCAN-3c energy "
                    "(within the already GFN2-xTB-filtered set).")
    parser.add_argument("--population-table", required=True,
                         help="Layer A's population_table.csv (from "
                              "filter_boltzmann_population.py)")
    parser.add_argument("--conformer-jobs", required=True,
                         help="Folder of job_NNNN/ dirs (from fanout_conformer_jobs.py), "
                              "job_0001 corresponds to the 1st surviving conformer "
                              "in population_filtered.xyz, in order")
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument("--temp", type=float, default=298.15)
    args = parser.parse_args()

    with open(args.population_table) as f:
        rows = [r for r in csv.DictReader(f) if r["kept"] == "1"]

    job_dirs = sorted(
        d for d in os.listdir(args.conformer_jobs)
        if d.startswith("job_")
    )
    if len(job_dirs) != len(rows):
        raise SystemExit(
            f"Mismatch: {len(rows)} kept conformers in population table "
            f"vs {len(job_dirs)} job directories -- these must come from "
            f"the same Layer A run.")

    energies_hartree = []
    missing = []
    for row, jd in zip(rows, job_dirs):
        e = read_last_energy(os.path.join(args.conformer_jobs, jd))
        if e is None:
            missing.append(jd)
        energies_hartree.append(e)

    if missing:
        raise SystemExit(
            f"{len(missing)} job(s) have no readable 'energy' file: "
            f"{missing[:5]}{'...' if len(missing) > 5 else ''} -- "
            f"fix or exclude these before re-ranking.")

    degeneracies = [int(row["degeneracy"]) for row in rows]

    e_min = min(energies_hartree)
    rel_kcal = [(e - e_min) * HARTREE_TO_KCAL for e in energies_hartree]

    kt = KB_KCAL_PER_MOL_K * args.temp
    unnorm = [deg * pow(2.718281828459045, -e / kt)
              for e, deg in zip(rel_kcal, degeneracies)]
    total = sum(unnorm)
    weights = [u / total for u in unnorm]

    combined = sorted(
        zip(rows, job_dirs, energies_hartree, rel_kcal, weights),
        key=lambda t: t[3],
    )

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "crest_index", "job_dir", "gfn2_relative_energy_kcal",
            "gfn2_population", "r2scan3c_energy_hartree",
            "r2scan3c_relative_energy_kcal", "r2scan3c_population",
        ])
        for row, jd, e, rel, w in combined:
            writer.writerow([
                row["crest_index"], jd, row["relative_energy_kcal"],
                row["population"], f"{e:.8f}", f"{rel:.4f}", f"{w:.6f}",
            ])

    top = combined[0]
    print(f"Re-ranked {len(combined)} Layer A survivors by r2SCAN-3c energy "
          f"(T={args.temp} K).")
    print(f"  new top conformer: {top[1]} (was crest_index={top[0]['crest_index']}, "
          f"GFN2-xTB relative E={top[0]['relative_energy_kcal']} kcal/mol) "
          f"-> r2SCAN-3c population {top[4]:.1%}")
    old_top_job = job_dirs[0]
    if top[1] != old_top_job:
        print(f"  NOTE: ranking changed -- GFN2-xTB's top conformer was "
              f"{old_top_job}, r2SCAN-3c's top conformer is {top[1]}")
    else:
        print("  GFN2-xTB and r2SCAN-3c agree on the top conformer.")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
