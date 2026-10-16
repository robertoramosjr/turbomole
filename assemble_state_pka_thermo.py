#!/usr/bin/env python3
"""
assemble_state_pka_thermo.py

Assembles the population-weighted G_gas and G_solv for one protonation
state, over its labeled (real HSE06 + freeh + COSMO) subset --
population weights renormalized WITHIN that subset, since it only
covers a fraction of the full Layer A ensemble (protonation_and_conformer_investigation.md
Sec. 9.32-9.33).

Usage:
    python ~/work_turbomole/scripts/assemble_state_pka_thermo.py \
        --state-label neutral \
        --labeled-jobs layer_b/labeled_jobs \
        --freeh-csv layer_b/freeh_thermal_correction.csv \
        --cosmo-csv layer_b/cosmo_energy.csv \
        --population-table layer_a/population_table_dft.csv \
        --output layer_b/state_pka_thermo.dat
"""

import argparse
import csv
import os


def main():
    parser = argparse.ArgumentParser(
        description="Population-weighted G_gas/G_solv for one protonation "
                    "state's labeled subset.")
    parser.add_argument("--state-label", required=True)
    parser.add_argument("--labeled-jobs", required=True)
    parser.add_argument("--freeh-csv", required=True)
    parser.add_argument("--cosmo-csv", required=True)
    parser.add_argument("--population-table", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.freeh_csv) as f:
        freeh = {r["job_dir"]: float(r["g_thermal_correction_hartree"])
                 for r in csv.DictReader(f)}
    with open(args.cosmo_csv) as f:
        cosmo = {r["job_dir"]: float(r["e_cosmo_oc_hartree"])
                 for r in csv.DictReader(f)}
    with open(args.population_table) as f:
        pop = {r["job_dir"]: float(r["r2scan3c_population"])
               for r in csv.DictReader(f)}

    common = sorted(set(freeh) & set(cosmo) & set(pop))
    HARTREE_TO_KCAL = 627.509474

    gas_energy = {}
    for jd in common:
        energy_path = os.path.join(args.labeled_jobs, jd, "energy")
        with open(energy_path) as f:
            lines = [l for l in f if l.strip() and not l.strip().startswith("$")]
        gas_energy[jd] = float(lines[-1].split()[1])

    w_sum = sum(pop[jd] for jd in common)
    g_gas_avg = sum(pop[jd] * (gas_energy[jd] + freeh[jd]) for jd in common) / w_sum
    g_solv_avg = sum(pop[jd] * (cosmo[jd] - gas_energy[jd]) for jd in common) / w_sum

    g_gas_kcal = g_gas_avg * HARTREE_TO_KCAL
    g_solv_kcal = g_solv_avg * HARTREE_TO_KCAL

    print(f"=== {args.state_label} (N={len(common)}, renormalized population "
          f"weight sum={w_sum:.4f}) ===")
    print(f"Population-weighted G_gas  = {g_gas_avg:.6f} Ha ({g_gas_kcal:.3f} kcal/mol)")
    print(f"Population-weighted G_solv = {g_solv_avg:.6f} Ha ({g_solv_kcal:.3f} kcal/mol)")

    with open(args.output, "w") as f:
        f.write("descriptor state,n,g_gas_hartree,g_gas_kcal,g_solv_hartree,g_solv_kcal\n")
        f.write(f"{args.state_label} {len(common)} {g_gas_avg:.6f} {g_gas_kcal:.3f} "
                f"{g_solv_avg:.6f} {g_solv_kcal:.3f}\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
