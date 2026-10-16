#!/usr/bin/env python3
"""
compute_pka_thermodynamic_cycle.py

Assembles the aqueous deprotonation free energy and predicted pKa for
an acid/conjugate-base pair, from HSE06/def2-TZVP gas-phase electronic
energies + freeh thermal corrections + COSMO(water) single-point
energies (protonation_and_conformer_investigation.md Sec. 9.32-9.33).

Thermodynamic cycle:
    G_gas(species)  = E_elec_gas(species) + G_thermal_correction(species)
    G_solv(species) = E_elec_COSMO(species) - E_elec_gas(species)
                      (single-point COSMO-vs-gas estimate of the
                      solvation free energy -- the thermal correction
                      itself is not recomputed in solution, a standard
                      simplification for this class of cheap implicit-
                      solvent cycle)
    deltaG_aq = [G_gas(base) + G_solv(base)]
              + G_proton_aq
              - [G_gas(acid) + G_solv(acid)]

G_proton_aq is the literature absolute free energy of the aqueous
proton (1 atm gas -> 1 M aqueous, combined so no separate standard-
state correction is needed for H+ itself): G(H+,gas,298K,1atm) =
-6.28 kcal/mol (ideal-gas translational free energy, Sackur-Tetrode --
this is the same value used throughout the computational pKa
literature, e.g. Kelly, Cramer & Truhlar) + deltaG_solv(H+) = -265.9
kcal/mol (Tissandier et al., J. Phys. Chem. A 1998, 102, 7787 --
cluster-pair approximation, the consensus literature value) =
-272.18 kcal/mol total, entered here as a single combined constant.

pKa_predicted = deltaG_aq / (R*T*ln(10))

Reports the RESIDUAL offset needed to match a given experimental pKa
-- this is the actual diagnostic (Sec. 9.32 Task 3/4): a residual
consistent between two different acid/base pairs at the same level of
theory is evidence of a systematic, transferable error (missing
outlying-charge effects, cavity parameterization, DFT error, etc.);
a residual that differs a lot between systems means the correction
does NOT transfer and should not be applied blindly elsewhere.

Usage:
    python ~/work_turbomole/scripts/compute_pka_thermodynamic_cycle.py \
        --label "artepilin_c_carboxyl" \
        --acid-e-gas -1234.567890 --acid-e-cosmo -1234.598765 \
        --acid-g-thermal-hartree 0.343116 \
        --base-e-gas -1234.111111 --base-e-cosmo -1234.155555 \
        --base-g-thermal-hartree 0.328920 \
        --experimental-pka 4.65 \
        --output pka_cycle_carboxyl.dat
"""

import argparse
import math

HARTREE_TO_KCAL = 627.509474
R_KCAL_MOL_K = 1.987204e-3
T_KELVIN = 298.15
G_PROTON_GAS_KCAL = -6.28       # ideal-gas translational free energy, Sackur-Tetrode
G_SOLV_PROTON_KCAL = -265.9     # Tissandier et al. 1998, cluster-pair approximation
G_PROTON_AQ_KCAL = G_PROTON_GAS_KCAL + G_SOLV_PROTON_KCAL


def compute_cycle(acid_e_gas, acid_e_cosmo, acid_g_thermal, base_e_gas,
                  base_e_cosmo, base_g_thermal):
    """Return (g_gas_acid, g_gas_base, g_solv_acid, g_solv_base,
    delta_g_aq, pka_predicted); energies in Hartree in, kcal/mol out."""
    g_gas_acid = (acid_e_gas + acid_g_thermal) * HARTREE_TO_KCAL
    g_gas_base = (base_e_gas + base_g_thermal) * HARTREE_TO_KCAL
    g_solv_acid = (acid_e_cosmo - acid_e_gas) * HARTREE_TO_KCAL
    g_solv_base = (base_e_cosmo - base_e_gas) * HARTREE_TO_KCAL
    delta_g_aq = (g_gas_base + g_solv_base) + G_PROTON_AQ_KCAL - (g_gas_acid + g_solv_acid)
    pka = delta_g_aq / (R_KCAL_MOL_K * T_KELVIN * math.log(10))
    return g_gas_acid, g_gas_base, g_solv_acid, g_solv_base, delta_g_aq, pka


def main():
    parser = argparse.ArgumentParser(
        description="Compute the aqueous deprotonation free energy and "
                    "predicted pKa from gas-phase + freeh + COSMO energies, "
                    "and the residual offset against a known experimental pKa.")
    parser.add_argument("--label", required=True)
    parser.add_argument("--acid-e-gas", type=float, required=True, help="Hartree")
    parser.add_argument("--acid-e-cosmo", type=float, required=True, help="Hartree, with OC correction")
    parser.add_argument("--acid-g-thermal-hartree", type=float, required=True)
    parser.add_argument("--base-e-gas", type=float, required=True, help="Hartree")
    parser.add_argument("--base-e-cosmo", type=float, required=True, help="Hartree, with OC correction")
    parser.add_argument("--base-g-thermal-hartree", type=float, required=True)
    parser.add_argument("--experimental-pka", type=float, default=None)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    g_gas_acid, g_gas_base, g_solv_acid, g_solv_base, delta_g_aq, pka_predicted = compute_cycle(
        args.acid_e_gas, args.acid_e_cosmo, args.acid_g_thermal_hartree,
        args.base_e_gas, args.base_e_cosmo, args.base_g_thermal_hartree)

    print(f"=== {args.label} ===")
    print(f"G_gas(acid)  = {g_gas_acid:.3f} kcal/mol")
    print(f"G_gas(base)  = {g_gas_base:.3f} kcal/mol")
    print(f"G_solv(acid) = {g_solv_acid:.3f} kcal/mol")
    print(f"G_solv(base) = {g_solv_base:.3f} kcal/mol")
    print(f"G_proton(aq) = {G_PROTON_AQ_KCAL:.3f} kcal/mol "
          f"(gas {G_PROTON_GAS_KCAL:.2f} + solv {G_SOLV_PROTON_KCAL:.2f}, "
          f"Tissandier et al. 1998)")
    print(f"\ndeltaG_aq(deprotonation) = {delta_g_aq:.3f} kcal/mol")
    print(f"pKa_predicted = {pka_predicted:.3f}")

    lines = [
        f"descriptor label,g_gas_acid_kcal,g_gas_base_kcal,g_solv_acid_kcal,"
        f"g_solv_base_kcal,g_proton_aq_kcal,delta_g_aq_kcal,pka_predicted",
    ]
    row = (f"{args.label} {g_gas_acid:.3f} {g_gas_base:.3f} {g_solv_acid:.3f} "
           f"{g_solv_base:.3f} {G_PROTON_AQ_KCAL:.3f} {delta_g_aq:.3f} "
           f"{pka_predicted:.3f}")

    if args.experimental_pka is not None:
        pka_residual = args.experimental_pka - pka_predicted
        g_residual_kcal = pka_residual * R_KCAL_MOL_K * T_KELVIN * math.log(10)
        print(f"\nExperimental pKa = {args.experimental_pka:.3f}")
        print(f"Residual offset needed to match: {pka_residual:+.3f} pKa units "
              f"({g_residual_kcal:+.3f} kcal/mol)")
        row += f",{args.experimental_pka:.3f},{pka_residual:.3f},{g_residual_kcal:.3f}"
        lines[0] += ",experimental_pka,pka_residual,g_residual_kcal"

    with open(args.output, "w") as f:
        f.write(lines[0] + "\n" + row + "\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
