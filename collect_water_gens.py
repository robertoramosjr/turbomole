#!/usr/bin/env python3
"""Free-energy ensemble of the aqueous conformers: E_HSE06/COSMO (W4) + G_thermal (xtb --hess).

For every conformer with a W4 single point (w4_sp/<state>/<cNNNN>/W4_DONE) the script reads
  E_W4  : "Total energy + OC corr." of ridft.out (HSE06-D3(BJ)/def2-TZVP, COSMO water, with the
          outlying-charge correction, as in the production cycle; Eh)
  G_th  : "G(RRHO) contrib." of xtb_gtherm_w4/<state>/<cNNNN>/hess.out (GFN2-xTB mRRHO; Eh)
  nimag : "# imaginary freq." of the same file
  g_i   : CREST rotamer degeneracy, line cNNNN of water_ensemble/<state>/cre_members
  E_r2  : r2SCAN-3c/COSMO energy of the re-optimised geometry (reopt/<state>/<cNNNN>/energy)
and builds, per state, G_i = E_W4 + G_th, the weights g_i exp(-G_i/kT) and, for comparison,
g_i exp(-E_r2/kT) (the interim production weights).  Conformers missing G_th are reported and
left out of the G ensemble (their E_W4 stays in the table).

G_ens = G_min - kT ln sum_i g_i exp(-(G_i - G_min)/kT)   (kcal/mol, 298.15 K, relative to
nothing: absolute values in Eh are written to the CSV).  Raw (uncalibrated) pKa per pair uses
G(H+,aq) = -272.18 kcal/mol; the class offsets are NOT applied here.

Outputs: <out_prefix>_conformers.csv, <out_prefix>_states.csv (and a stdout summary).

Usage:
  python collect_water_gens.py --root ~/work_turbomole/artepillin_C_studies \
      --out-prefix paper_2/Figures/data/water_gens
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
PAIRS = [("N_to_A", "N", "A"), ("N_to_B", "N", "B"), ("A_to_D", "A", "D"), ("B_to_D", "B", "D")]
RE_E = re.compile(r"Total energy \+ OC corr\.\s*=\s*(-?\d+\.\d+)")
RE_G = re.compile(r"G\(RRHO\) contrib\.\s+(-?\d+\.\d+)\s+Eh")
RE_IMAG = re.compile(r"# imaginary freq\.\s+(\d+)")


def last_match(path, regex):
    if not os.path.isfile(path):
        return None
    m = regex.findall(open(path, errors="ignore").read())
    return m[-1] if m else None


def read_degeneracy(path):
    """Degeneracy of CREST structure k (1-based) = first column of line k+1 of cre_members."""
    lines = [l.split() for l in open(path) if l.strip()]
    return [int(l[0]) for l in lines if len(l) >= 3]


def read_energy_file(path):
    if not os.path.isfile(path):
        return None
    rows = [l.split() for l in open(path) if l.strip() and not l.startswith("$")]
    return float(rows[-1][1]) if rows else None


def collect(root, label):
    state = STATES[label]
    deg = read_degeneracy(os.path.join(root, "water_ensemble", state, "cre_members"))
    base = os.path.join(root, "w4_sp", state)
    recs = []
    for job in sorted(os.listdir(base)):
        d = os.path.join(base, job)
        if not os.path.isfile(os.path.join(d, "W4_DONE")):
            continue
        idx = int(job[1:])
        e_w4 = last_match(os.path.join(d, "ridft.out"), RE_E)
        xd = os.path.join(root, "xtb_gtherm_w4", state, job, "hess.out")
        g_th, nim = last_match(xd, RE_G), last_match(xd, RE_IMAG)
        recs.append(dict(state=label, job=job, crest=idx, deg=deg[idx - 1],
                         e_w4=float(e_w4), g_th=None if g_th is None else float(g_th),
                         nimag=None if nim is None else int(nim),
                         e_r2=read_energy_file(os.path.join(root, "reopt", state, job, "energy"))))
    return recs


def boltzmann(recs, key):
    """Normalised weights g_i exp(-key_i/kT) (key in Eh)."""
    ref = min(r[key] for r in recs)
    w = [r["deg"] * math.exp(-(r[key] - ref) * HARTREE_TO_KCAL / KT) for r in recs]
    z = sum(w)
    return [x / z for x in w]


def g_ens(recs):
    """Ensemble free energy (Eh): G_min - kT ln sum g_i exp(-(G_i - G_min)/kT)."""
    gmin = min(r["g"] for r in recs)
    z = sum(r["deg"] * math.exp(-(r["g"] - gmin) * HARTREE_TO_KCAL / KT) for r in recs)
    return gmin - KT * math.log(z) / HARTREE_TO_KCAL


def neff(w):
    return 1.0 / sum(x * x for x in w)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    out = os.path.abspath(os.path.expanduser(a.out_prefix))

    allrecs, summary, gens = {}, [], {}
    for label in STATES:
        recs = collect(root, label)
        miss = [r["job"] for r in recs if r["g_th"] is None]
        bad = [r["job"] for r in recs if r["nimag"] not in (0, None)]
        use = [r for r in recs if r["g_th"] is not None and r["nimag"] == 0]
        for r in use:
            r["g"] = r["e_w4"] + r["g_th"]
        wg = boltzmann(use, "g")
        we = boltzmann(use, "e_w4")
        w2 = boltzmann([r for r in use if r["e_r2"] is not None], "e_r2") if all(r["e_r2"] is not None for r in use) else None
        for r, x in zip(use, wg):
            r["w_G"] = x
        for r, x in zip(use, we):
            r["w_EW4"] = x
        if w2:
            for r, x in zip(use, w2):
                r["w_Er2"] = x
        tv = lambda p, q: 0.5 * sum(abs(u - v) for u, v in zip(p, q))
        gens[label] = g_ens(use)
        summary.append(dict(state=label, n_w4=len(recs), n_used=len(use), n_missing_xtb=len(miss),
                            n_imag=len(bad), neff_G=neff(wg), neff_EW4=neff(we),
                            neff_Er2=neff(w2) if w2 else "", tv_G_vs_EW4=tv(wg, we),
                            tv_G_vs_Er2=tv(wg, w2) if w2 else "", G_ens_Eh=gens[label],
                            missing=";".join(miss), imaginary=";".join(bad)))
        allrecs[label] = use + [r for r in recs if r not in use]

    os.makedirs(os.path.dirname(out), exist_ok=True)
    cols = ["state", "job", "crest", "deg", "e_w4", "g_th", "nimag", "e_r2", "g", "w_G", "w_EW4", "w_Er2"]
    with open(out + "_conformers.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for label in STATES:
            for r in allrecs[label]:
                w.writerow(r)
    with open(out + "_states.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]))
        w.writeheader()
        for s in summary:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in s.items()})
    for s in summary:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items() if k != "G_ens_Eh"})
    for name, acid, base in PAIRS:
        raw = ((gens[base] - gens[acid]) * HARTREE_TO_KCAL + G_PROTON_AQ) / (KT * math.log(10))
        print(f"{name}: raw pKa (G_ens, water, no class offset) = {raw:.3f}")
    print("wrote", out + "_conformers.csv", out + "_states.csv")


if __name__ == "__main__":
    main()
