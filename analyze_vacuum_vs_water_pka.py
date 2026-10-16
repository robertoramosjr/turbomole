#!/usr/bin/env python3
"""
analyze_vacuum_vs_water_pka.py

Decomposes the vacuum-vs-water pKa discrepancy of artepilin C (Sec. 9.39-9.40,
compare_vacuum_vs_water_pka.py). That script mixes definitions: variant A weights
by the r2SCAN-3c (vacuum, degeneracy-weighted) population, while B/C/D Boltzmann-
weight the aqueous free energy WITHOUT the CREST rotamer degeneracy, and A-C use
the population-averaged <G> instead of the ensemble free energy
G_ens = -kT ln sum_i g_i exp(-G_i/kT). This script recomputes the pKa with every
combination so each effect can be isolated:

  weights : pop   = r2SCAN-3c vacuum population (production)
            gaq   = Boltzmann on G_aq, no degeneracy   (old B/C)
            gaqg  = Boltzmann on G_aq WITH degeneracy g_i (correct statistical weight)
  state G : mean  = sum w_i G_i              (production style)
            ens   = -kT ln sum g_i exp(-G_i/kT) (includes rotamer degeneracy and mixing entropy;
                    for subsets, divided by the covered r2SCAN-3c population fraction)
  set     : old (labelled subset) | all (old + water-guided new conformers)

G_aq(i) = E_COSMO(OC)(i) + G_thermal(i) [E_gas cancels in the cycle, as in the
production cycle]. Output: per-conformer table and a pKa table (raw, i.e. before the
class offset, and calibrated with the default offsets).

Usage:
    python analyze_vacuum_vs_water_pka.py --project-root ~/work_turbomole/artepillin_C_studies \
        --outprefix vac_vs_water_decomposition
"""

import argparse
import csv
import math
import os
import re

HARTREE_TO_KCAL = 627.509474
KT = 0.0019872041 * 298.15
G_PROTON_AQ = -272.18
STATES = {"N": "neutral", "A": "monoanion_a", "B": "monoanion_b", "D": "dianion"}
EXT_DIR = {"A": "top30_ext", "B": "water_ext", "D": "water_ext"}
PAIRS = [("N_to_A", "N", "A", "0to-1"), ("N_to_B", "N", "B", "0to-1"),
         ("A_to_D", "A", "D", "-1to-2"), ("B_to_D", "B", "D", "-1to-2")]


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def gaq(layer_b):
    fr = {r["job_dir"]: float(r["g_thermal_correction_hartree"])
          for r in rows(os.path.join(layer_b, "freeh_thermal_correction.csv"))}
    co = {r["job_dir"]: float(r["e_cosmo_oc_hartree"]) for r in rows(os.path.join(layer_b, "cosmo_energy.csv"))}
    return {j: (co[j] + fr[j]) * HARTREE_TO_KCAL for j in set(fr) & set(co)}


def egas_hse06(job_dir):
    """HSE06 gas-phase total energy (Hartree) from the last 'total energy' line of ridft.out."""
    p = os.path.join(job_dir, "ridft.out")
    val = None
    if os.path.exists(p):
        with open(p, errors="ignore") as f:
            for ln in f:
                m = re.search(r"total energy\s*=\s*(-?\d+\.\d+)", ln)
                if m:
                    val = float(m.group(1))
    return val


def cosmo_hartree(layer_b):
    return {r["job_dir"]: float(r["e_cosmo_oc_hartree"]) for r in rows(os.path.join(layer_b, "cosmo_energy.csv"))}


def load_state(root, s):
    base = os.path.join(root, f"artepillin_c_{STATES[s]}_d3_hse06")
    with open(os.path.join(base, "cre_members")) as f:
        deg = [int(ln.split()[0]) for ln in f if ln.strip() and ln.split()[0].isdigit()]
    # first line of cre_members is the total count -> keep only per-conformer lines
    with open(os.path.join(base, "cre_members")) as f:
        lines = [ln.split() for ln in f if ln.strip()]
    deg = [int(l[0]) for l in lines if len(l) >= 3]
    pop_dft = {r["job_dir"]: r for r in rows(os.path.join(base, "layer_a", "population_table_dft.csv"))}
    pop_total = sum(float(r["r2scan3c_population"]) for r in pop_dft.values())
    old = gaq(os.path.join(base, "layer_b"))
    recs = []
    for j, g in old.items():
        if j in pop_dft:
            ci = int(pop_dft[j]["crest_index"])
            eg = egas_hse06(os.path.join(base, "layer_b", "labeled_jobs", j))
            ec = cosmo_hartree(os.path.join(base, "layer_b")).get(j)
            recs.append(dict(set="old", job=j, crest=ci, g=g, deg=deg[ci - 1], egas=eg, ecosmo=ec,
                             pop=float(pop_dft[j]["r2scan3c_population"]) / pop_total))
    if s in EXT_DIR:
        ed = os.path.join(base, "layer_b", EXT_DIR[s])
        if os.path.exists(os.path.join(ed, "freeh_thermal_correction.csv")):
            mp = {r["job_dir"]: int(r["crest_index"]) for r in rows(os.path.join(base, "layer_a", "water_ext", "map.csv"))} \
                if s != "A" else None
            if s == "A":
                # top30_ext: job names carry their own map file
                mpf = os.path.join(base, "layer_b", "top30_ext", "map.csv")
                mp = {r["job_dir"]: int(r["crest_index"]) for r in rows(mpf)} if os.path.exists(mpf) else {}
            ec_new = cosmo_hartree(ed)
            sub = "labeled_jobs"
            for j, g in gaq(ed).items():
                ci = mp.get(j) if s != "A" else int(j.split("_")[1])  # A: job number = CREST index
                eg = egas_hse06(os.path.join(ed, sub, j))
                recs.append(dict(set="new", job=j, crest=ci, g=g, deg=deg[ci - 1] if ci else 1, egas=eg,
                                 ecosmo=ec_new.get(j),
                                 pop=(float(pop_dft[j]["r2scan3c_population"]) / pop_total) if (s == "A" and j in pop_dft) else 0.0))
    return recs, pop_total


def state_G(recs, weights, kind, coverage=1.0):
    """weights in {'pop','gaq','gaqg'}; kind in {'mean','ens'}."""
    gmin = min(r["g"] for r in recs)
    if kind == "ens":
        z = sum(r["deg"] * math.exp(-(r["g"] - gmin) / KT) for r in recs)
        return gmin - KT * math.log(z / coverage)
    if weights == "pop":
        w = [r["pop"] for r in recs]
    elif weights == "gaq":
        w = [math.exp(-(r["g"] - gmin) / KT) for r in recs]
    else:
        w = [r["deg"] * math.exp(-(r["g"] - gmin) / KT) for r in recs]
    tot = sum(w)
    return sum(wi * r["g"] for wi, r in zip(w, recs)) / tot


def pka(ga, gb):
    return (gb + G_PROTON_AQ - ga) / (KT * math.log(10))


def main():
    ap = argparse.ArgumentParser(description="Decompose the vacuum-vs-water pKa discrepancy.")
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--outprefix", required=True)
    ap.add_argument("--offset-0to-1", type=float, default=-3.785)
    ap.add_argument("--offset-1to-2", type=float, default=-3.727)
    a = ap.parse_args()
    root = os.path.expanduser(a.project_root)
    offs = {"0to-1": a.offset_0to_1, "-1to-2": a.offset_1to_2}
    data = {s: load_state(root, s) for s in STATES}

    # per-conformer table (relative G_aq within each state)
    with open(a.outprefix + "_conformers.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["state", "set", "job", "crest_index", "degeneracy", "G_aq_rel_kcal", "vac_pop_renorm",
                    "dE_gas_HSE06_kcal", "dG_solv_COSMO_kcal"])
        for s, (recs, _) in data.items():
            gmin = min(r["g"] for r in recs)
            ref = min(recs, key=lambda r: r["g"])
            for r in sorted(recs, key=lambda r: r["g"]):
                dg = dsol = ""
                if r.get("egas") is not None and ref.get("egas") is not None and r.get("ecosmo") is not None:
                    dg = f"{(r['egas'] - ref['egas']) * HARTREE_TO_KCAL:.3f}"
                    dsol = f"{((r['ecosmo'] - r['egas']) - (ref['ecosmo'] - ref['egas'])) * HARTREE_TO_KCAL:.3f}"
                w.writerow([s, r["set"], r["job"], r["crest"], r["deg"], f"{r['g'] - gmin:.3f}", f"{r['pop']:.4f}", dg, dsol])

    schemes = [("pop", "mean", "old"), ("gaq", "mean", "old"), ("gaqg", "mean", "old"), ("gaqg", "ens", "old"),
               ("gaq", "mean", "all"), ("gaqg", "mean", "all"), ("gaqg", "ens", "all")]
    G = {}
    for s, (recs, ptot) in data.items():
        old = [r for r in recs if r["set"] == "old"]
        cov_old = sum(r["pop"] for r in old)  # fraction of the vacuum-95% population that has a G_aq
        for sc in schemes:
            sub = old if sc[2] == "old" else recs
            if sc[0] == "pop" and sc[2] == "all":
                continue
            # ens over a subset: divide by covered fraction of the r2SCAN-3c population (best available estimate)
            # coverage of the vacuum-95% population by the conformers that carry a G_aq (only meaningful for N, A;
            # B and D are fully covered by construction, the water-guided extras lie outside that population)
            cov = (sum(r["pop"] for r in sub) if s in "NA" else 1.0) if sc[1] == "ens" else 1.0
            G[(s, sc)] = state_G(sub, sc[0], sc[1], cov)
        print(f"{s}: n_old={len(old)} n_new={len(recs) - len(old)} vac-pop coverage of old subset={cov_old:.3f}")
    with open(a.outprefix + "_pka.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pair", "weights", "state_G", "set", "raw_pKa", "calib_pKa", "calib_minus_production"])
        for name, acid, base, cls in PAIRS:
            prod = pka(G[(acid, schemes[0])], G[(base, schemes[0])]) + offs[cls]
            for sc in schemes:
                if (acid, sc) in G and (base, sc) in G:
                    raw = pka(G[(acid, sc)], G[(base, sc)])
                    w.writerow([name, sc[0], sc[1], sc[2], f"{raw:.3f}", f"{raw + offs[cls]:.3f}",
                                f"{raw + offs[cls] - prod:+.3f}"])
                    print(f"{name:7s} w={sc[0]:5s} G={sc[1]:4s} set={sc[2]:3s} raw={raw:7.3f} calib={raw + offs[cls]:7.3f} "
                          f"d_vs_production={raw + offs[cls] - prod:+.3f}")
    print("wrote", a.outprefix + "_conformers.csv", a.outprefix + "_pka.csv")


if __name__ == "__main__":
    main()
