#!/usr/bin/env python3
"""
calibrate_pka_pair.py

One-shot calibration of a deprotonation pair
(protonation_and_conformer_investigation.md Sec. 9.33/9.34): computes
the raw thermodynamic-cycle pKa of a small REFERENCE compound
(acid/base pair with a well-established aqueous pKa), takes the
residual against its experimental pKa as a systematic offset, and
applies it to the artepilin C TARGET pair, whose per-state
population-weighted G_gas/G_solv come from assemble_state_pka_thermo.py.

    pKa_target_calibrated = pKa_target_raw + residual_reference,
    residual_reference    = pKa_exp(reference) - pKa_raw(reference)

The cycle itself is compute_pka_thermodynamic_cycle.compute_cycle (one
formula for everything -- G_proton_aq = -272.18 kcal/mol, Tissandier
1998 + Sackur-Tetrode). Caveats printed with every result:
  * N=1 reference per chemical class: the offset's transferability to a
    DIFFERENT class (carboxyl -> phenol) is untested unless a second
    reference of the target class exists.
  * --target-sconf-shift-kcal adds the difference in ensemble
    conformational free energy (-T*S_conf, from the FULL Layer A
    populations incl. CREST rotamer degeneracy) between the base and
    acid states as a sensitivity term; default 0 keeps the plain
    population-weighted-average definition used so far.

Usage:
    python ~/work_turbomole/scripts/calibrate_pka_pair.py \
        --label phenol_pair \
        --ref-acid-dir pka_reference/phenol/acid \
        --ref-base-dir pka_reference/phenol/base \
        --ref-exp-pka 9.99 \
        --target-acid-thermo artepillin_c_neutral_d3_hse06/layer_b/state_pka_thermo.dat \
        --target-base-thermo artepillin_c_monoanion_b_d3_hse06/layer_b/state_pka_thermo.dat \
        --target-exp-pka 4.65 --target-sconf-shift-kcal 0.303 \
        --output pka_calibration_phenol_pair.dat
"""

import argparse
import csv
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compute_pka_thermodynamic_cycle import (  # noqa: E402
    compute_cycle, R_KCAL_MOL_K, T_KELVIN,
)

OC_RE = re.compile(r"Total energy \+ OC corr\.\s*=\s*([-\d.]+)")


def read_reference_species(species_dir):
    with open(os.path.join(species_dir, "energy")) as f:
        rows = [l for l in f if l.strip() and not l.strip().startswith("$")]
    e_gas = float(rows[-1].split()[1])
    with open(os.path.join(species_dir + "_cosmo", "ridft_cosmo.out")) as f:
        m = OC_RE.search(f.read())
    if not m:
        raise SystemExit(f"no OC-corrected COSMO energy in {species_dir}_cosmo")
    e_cosmo = float(m.group(1))
    with open(os.path.join(species_dir, "freeh_thermal_correction.csv")) as f:
        g_th = float(next(csv.DictReader(f))["g_thermal_correction_hartree"])
    return e_gas, e_cosmo, g_th


def read_state_thermo(path):
    with open(path) as f:
        lines = [l.split() for l in f if l.strip() and not l.startswith("descriptor")]
    state, n, g_gas_ha, g_gas_kcal, g_solv_ha, g_solv_kcal = lines[0]
    return float(g_gas_ha), float(g_solv_ha), int(n)


def main():
    parser = argparse.ArgumentParser(
        description="Calibrate a deprotonation pair against a reference "
                    "compound's residual and apply it to the artepilin C pair.")
    parser.add_argument("--label", required=True)
    parser.add_argument("--ref-acid-dir", required=True)
    parser.add_argument("--ref-base-dir", required=True)
    parser.add_argument("--ref-exp-pka", type=float, required=True)
    parser.add_argument("--target-acid-thermo", required=True)
    parser.add_argument("--target-base-thermo", required=True)
    parser.add_argument("--target-exp-pka", type=float, default=None,
                         help="Experimental pKa of the target pair, if one exists "
                              "(only for reporting the final residual)")
    parser.add_argument("--target-sconf-shift-kcal", type=float, default=0.0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    ra = read_reference_species(args.ref_acid_dir)
    rb = read_reference_species(args.ref_base_dir)
    *_, ref_dg, ref_pka_raw = compute_cycle(ra[0], ra[1], ra[2], rb[0], rb[1], rb[2])
    ref_res = args.ref_exp_pka - ref_pka_raw
    print(f"=== reference ({os.path.basename(os.path.dirname(args.ref_acid_dir.rstrip('/')))}) ===")
    print(f"raw pKa = {ref_pka_raw:.3f}, experimental = {args.ref_exp_pka:.3f}, "
          f"residual = {ref_res:+.3f} pKa ({ref_res * R_KCAL_MOL_K * T_KELVIN * math.log(10):+.3f} kcal/mol)")

    ta_g, ta_s, na = read_state_thermo(args.target_acid_thermo)
    tb_g, tb_s, nb = read_state_thermo(args.target_base_thermo)
    *_, tg_dg, tg_pka_raw = compute_cycle(ta_g, ta_g + ta_s, 0.0, tb_g, tb_g + tb_s, 0.0)
    shift_pka = args.target_sconf_shift_kcal / (R_KCAL_MOL_K * T_KELVIN * math.log(10))
    tg_pka_raw_s = tg_pka_raw + shift_pka
    cal = tg_pka_raw + ref_res
    cal_s = tg_pka_raw_s + ref_res
    print(f"\n=== target pair (N_acid={na}, N_base={nb} labeled conformers) ===")
    print(f"raw pKa = {tg_pka_raw:.3f}  ->  calibrated = {cal:.3f}")
    if args.target_sconf_shift_kcal:
        print(f"with S_conf sensitivity ({args.target_sconf_shift_kcal:+.3f} kcal/mol = "
              f"{shift_pka:+.3f} pKa): raw {tg_pka_raw_s:.3f}  ->  calibrated {cal_s:.3f}")
    if args.target_exp_pka is not None:
        print(f"experimental = {args.target_exp_pka:.3f}; "
              f"calibrated - experimental = {cal - args.target_exp_pka:+.3f} pKa")
    print("\nCAVEAT: offset from ONE reference of this chemical class; "
          "transferability to the target is not independently tested.")

    with open(args.output, "w") as f:
        f.write("descriptor label,ref_pka_raw,ref_pka_exp,ref_residual,"
                "target_pka_raw,target_pka_calibrated,sconf_shift_pka,"
                "target_pka_calibrated_with_sconf\n")
        f.write(f"{args.label} {ref_pka_raw:.3f} {args.ref_exp_pka:.3f} {ref_res:.3f} "
                f"{tg_pka_raw:.3f} {cal:.3f} {shift_pka:.3f} {cal_s:.3f}\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
