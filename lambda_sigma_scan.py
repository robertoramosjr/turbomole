#!/usr/bin/env python3
"""Sensitivity of the aqueous-ensemble band maximum to the Gaussian broadening sigma.

Reads water_tddft/<state>/<cNNNN>/exspectrum (COSMO TD-DFT, length-gauge f) and the Boltzmann weights
w_G of water_gens_conformers.csv; builds eps(E) = 2.174e8/ (sigma*sqrt(2 pi)) ... in eV units, normalised
over the conformers that have an escf spectrum, and reports lambda_max (nm) for each sigma (eps scale differs from water_uv_summary.csv by a constant 1.51 and is not reported).
The sigma=0.4 eV row must reproduce water_uv_summary.csv (checked with --check).

Usage:
  python lambda_sigma_scan.py --root ~/work_turbomole/artepillin_C_studies \
      --sigmas 0.1 0.2 0.3 0.4 --out todo_v12/analysis/lambda_sigma_scan.csv [--check]
"""
import argparse
import csv
import os
import re

import numpy as np

STATES = {"N": "neutral", "A": "monoanion_a", "B": "monoanion_b", "D": "dianion"}
HC = 1239.84198  # eV nm
# eps [L mol-1 cm-1] = 2.174e8 * f / (sigma_cm-1 * sqrt(2 pi)) for a Gaussian of std sigma_cm-1
EV_TO_CM = 8065.5439


def read_exspectrum(path):
    out = []
    for line in open(path):
        m = re.match(r"\s*\d+\s+a\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)", line)
        if m:
            out.append((float(m.group(1)), float(m.group(3))))  # eV, f(length)
    return np.array(out)


def band(spectra, weights, sigma_ev, lam_lo=250.0, lam_hi=500.0):
    e = np.linspace(HC / lam_hi, HC / lam_lo, 4000)
    eps = np.zeros_like(e)
    tot = sum(weights)
    s_cm = sigma_ev * EV_TO_CM
    for sp, w in zip(spectra, weights):
        for ev, f in sp:
            eps += (w / tot) * 2.174e8 * f / (s_cm * np.sqrt(2 * np.pi)) * np.exp(-((e - ev) ** 2) / (2 * sigma_ev ** 2))
    i = int(np.argmax(eps))
    return HC / e[i], eps[i]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--sigmas", type=float, nargs="+", default=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4])
    ap.add_argument("--out", required=True)
    ap.add_argument("--check", action="store_true", help="print the sigma=0.4 eV values next to water_uv_summary.csv")
    a = ap.parse_args()
    root = os.path.expanduser(a.root)
    gens = list(csv.DictReader(open(os.path.join(root, "paper_2/Figures/data/water_gens_conformers.csv"))))
    rows = []
    for st, d in STATES.items():
        sp, w = [], []
        for r in gens:
            if r["state"] != st:
                continue
            p = os.path.join(root, "water_tddft", d, r["job"], "exspectrum")
            if os.path.isfile(p):
                sp.append(read_exspectrum(p))
                w.append(float(r["w_G"]))
        for s in a.sigmas:
            lam, eps = band(sp, w, s)
            rows.append({"state": st, "n_conf": len(sp), "sigma_eV": s, "lambda_max_nm": round(float(lam), 2)})
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        print(r)
    if a.check:
        print(open(os.path.join(root, "paper_2/Figures/data/water_uv_summary.csv")).read())


if __name__ == "__main__":
    main()
