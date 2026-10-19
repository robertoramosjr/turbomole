#!/usr/bin/env python3
"""Final pKa estimate for artepillin C with explicit uncertainty (Paper 2, sec:pka).

Inputs already on disk, nothing is computed here:
  * production microscopic values (gaqg weights, G_ens, water-extended set) from
    Figures/data/vac_vs_water_decomp_pka.csv;
  * class scatter of the calibration offsets (N=3 for 0->-1, N=2 for -1->-2) from
    pka_refs.csv and the raw reference cycle values;
  * an independent analogue estimate of the phenol step of the monoanion
    (2,6-dimethylphenol + [4-coumarate - phenol]), built only from reference
    pKa values already in pka_refs.csv.
Macroscopic constants follow from Ka1 = KNA + KNB and Ka1*Ka2 = KNA*KAD.
"""
import argparse
import csv
import math
import os

SOLVENT_U = 0.4  # model-of-solvent uncertainty, midpoint of the 0.3-0.5 range of the log
CLASS_SD = {"N_to_A": 0.406, "N_to_B": 0.406, "A_to_D": 1.152, "B_to_D": 1.152}


def production(path):
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            if (r["weights"], r["state_G"], r["set"]) == ("gaqg", "ens", "all"):
                out[r["pair"]] = float(r["calib_pKa"])
    return out


def macro(p):
    ka = lambda k: 10 ** -p[k]
    ka1 = ka("N_to_A") + ka("N_to_B")
    pka1 = -math.log10(ka1)
    pka2 = p["N_to_A"] + p["A_to_D"] - pka1
    return pka1, pka2


def quad(*u):
    return math.sqrt(sum(x * x for x in u))


def analogue(refs, coumarate_exp=None):
    e = {r["name"]: float(r["exp_pka"]) for r in refs}
    cou = e["4-coumarate_phenol"] if coumarate_exp is None else coumarate_exp
    return e["26-dimethylphenol"] + (cou - e["phenol"])


def run(root, out):
    p = production(os.path.join(root, "paper_2/Figures/data/vac_vs_water_decomp_pka.csv"))
    with open(os.path.join(root, "pka_refs.csv")) as f:
        refs = list(csv.DictReader(f))
    pka1, pka2 = macro(p)
    rows = [(k, p[k], quad(CLASS_SD[k], SOLVENT_U)) for k in ("N_to_A", "N_to_B", "A_to_D", "B_to_D")]
    rows += [("macro_pKa1", pka1, quad(CLASS_SD["N_to_A"], SOLVENT_U)),
             ("macro_pKa2", pka2, quad(CLASS_SD["A_to_D"], SOLVENT_U))]
    lo, hi = analogue(refs), analogue(refs, 9.92)
    rows.append(("analogue_pKa2", 0.5 * (lo + hi), 0.5))
    with open(out, "w") as f:
        f.write("quantity,pKa,uncertainty\n")
        for k, v, u in rows:
            f.write(f"{k},{v:.2f},{u:.2f}\n")
    return rows, (lo, hi)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", default=os.path.expanduser("~/work_turbomole/artepillin_C_studies"))
    ap.add_argument("--out", default="Figures/data/pka_final_estimate.csv")
    a = ap.parse_args()
    rows, rng = run(a.root, a.out)
    for r in rows:
        print("%-14s %6.2f +/- %.2f" % r)
    print("analogue range (coumarate 9.45 / 9.92): %.2f / %.2f" % rng)
