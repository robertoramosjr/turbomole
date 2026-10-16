#!/usr/bin/env python3
"""
compare_vacuum_vs_water_pka.py

Compares the artepilin C pKa obtained with the vacuum-CREST ensemble
(current production: labeled subset weighted by r2SCAN-3c population)
against ensembles that include the water-guided extra conformers
(ALPB rerank, Sec. 9.39-9.40), weighted by the aqueous free energy
G_aq(i) = E_COSMO(OC)(i) + G_thermal(i)   [E_gas cancels in the cycle].

Variants per protonation state (neutral has no water extension, so it
only differs by weighting):
  A  old labeled subset, r2SCAN-3c population   (= production, reproduces raw pKa)
  B  old labeled subset, Boltzmann on G_aq      (effect of weighting only)
  C  old + new conformers, Boltzmann on G_aq    (effect of ensemble extension)
  D  as C but ensemble free energy -kT ln sum exp(-G/kT) (adds conformational
     entropy; sensitivity, same spirit as --apply-sconf)
Population-weighted averages of G use arithmetic mean for A-C (as
assemble_state_pka_thermo.py does).

Usage:
    python compare_vacuum_vs_water_pka.py --project-root ~/work_turbomole/artepillin_C_studies \
        --output pka_vacuum_vs_water.csv
"""

import argparse
import csv
import math
import os

HARTREE_TO_KCAL = 627.509474
R = 0.0019872041
T = 298.15
G_PROTON_AQ = -272.18
KT = R * T

STATES = {
    "neutral": "artepillin_c_neutral_d3_hse06",
    "A": "artepillin_c_monoanion_a_d3_hse06",
    "B": "artepillin_c_monoanion_b_d3_hse06",
    "D": "artepillin_c_dianion_d3_hse06",
}
# subdirectory of layer_b holding the water-guided extra conformers
EXT_DIR = {"A": "top30_ext", "B": "water_ext", "D": "water_ext"}
PAIRS = [("N_to_A", "neutral", "A", "0to-1"), ("N_to_B", "neutral", "B", "0to-1"),
         ("A_to_D", "A", "D", "-1to-2"), ("B_to_D", "B", "D", "-1to-2")]


def read_col(path, col):
    with open(path) as f:
        return {r["job_dir"]: float(r[col]) for r in csv.DictReader(f)}


def g_aq_table(layer_b):
    freeh = read_col(os.path.join(layer_b, "freeh_thermal_correction.csv"),
                     "g_thermal_correction_hartree")
    cosmo = read_col(os.path.join(layer_b, "cosmo_energy.csv"), "e_cosmo_oc_hartree")
    return {j: (cosmo[j] + freeh[j]) * HARTREE_TO_KCAL for j in set(freeh) & set(cosmo)}


def boltzmann(g):
    gmin = min(g.values())
    w = {k: math.exp(-(v - gmin) / KT) for k, v in g.items()}
    z = sum(w.values())
    return {k: x / z for k, x in w.items()}


def state_free_energies(root, state):
    """Return dict variant -> G_aq (kcal/mol) for one state, plus sizes."""
    base = os.path.join(root, STATES[state])
    old = g_aq_table(os.path.join(base, "layer_b"))
    with open(os.path.join(base, "layer_a", "population_table_dft.csv")) as f:
        pop = {r["job_dir"]: float(r["r2scan3c_population"]) for r in csv.DictReader(f)}
    common = sorted(set(old) & set(pop))
    wsum = sum(pop[j] for j in common)
    out = {"A": sum(pop[j] * old[j] for j in common) / wsum}
    wb = boltzmann({j: old[j] for j in common})
    out["B"] = sum(wb[j] * old[j] for j in common)
    n_old, n_new = len(common), 0
    union = {("old", j): old[j] for j in common}
    if state in EXT_DIR:
        ext_dir = os.path.join(base, "layer_b", EXT_DIR[state])
        if os.path.exists(os.path.join(ext_dir, "freeh_thermal_correction.csv")):
            new = g_aq_table(ext_dir)
        else:
            print(f"WARNING: no freeh/cosmo csv in {ext_dir}; state {state} "
                  f"uses the old subset only (variant C == B)")
            new = {}
        n_new = len(new)
        union.update({("new", j): g for j, g in new.items()})
    wu = boltzmann(union)
    out["C"] = sum(wu[k] * union[k] for k in union)
    gmin = min(union.values())
    out["D"] = gmin - KT * math.log(sum(math.exp(-(g - gmin) / KT) for g in union.values()))
    return out, n_old, n_new


def pka(g_acid, g_base):
    return (g_base + G_PROTON_AQ - g_acid) / (KT * math.log(10))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--project-root", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--offset-0to-1", type=float, default=-3.785,
                   help="class offset (calibrated - raw), default from 5-ref summary")
    p.add_argument("--offset-1to-2", type=float, default=-3.727)
    a = p.parse_args()
    root = os.path.expanduser(a.project_root)
    offs = {"0to-1": a.offset_0to_1,
            "-1to-2": a.offset_1to_2}
    res, sizes = {}, {}
    for s in STATES:
        res[s], n_old, n_new = state_free_energies(root, s)
        sizes[s] = (n_old, n_new)
        print(f"state {s}: N_old={n_old} N_new={n_new}")
    rows = []
    for name, acid, base, cls in PAIRS:
        raw = {v: pka(res[acid][v], res[base][v]) for v in "ABCD"}
        row = {"pair": name}
        for v in "ABCD":
            row[f"raw_{v}"] = f"{raw[v]:.3f}"
            row[f"calib_{v}"] = f"{raw[v] + offs[cls]:.3f}"
        row["delta_C_minus_A"] = f"{raw['C'] - raw['A']:+.3f}"
        row["delta_B_minus_A"] = f"{raw['B'] - raw['A']:+.3f}"
        row["delta_D_minus_C"] = f"{raw['D'] - raw['C']:+.3f}"
        rows.append(row)
        print(f"{name}: raw A={raw['A']:.3f} B={raw['B']:.3f} C={raw['C']:.3f} D={raw['D']:.3f} | "
              f"calib A={raw['A'] + offs[cls]:.3f} C={raw['C'] + offs[cls]:.3f} | "
              f"dB={row['delta_B_minus_A']} dC={row['delta_C_minus_A']}")
    with open(a.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {a.output}")


if __name__ == "__main__":
    main()
