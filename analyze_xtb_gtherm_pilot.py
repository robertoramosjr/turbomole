#!/usr/bin/env python3
"""Compare xtb (GFN2, mRRHO) G_thermal with HSE06 freeh G_thermal, conformer by conformer.

Reads <root>/xtb_gtherm_pilot/<state>/<job>/{bhess,hess}.out ("G(RRHO) contrib." in Eh) and the HSE06
layer_b/freeh_thermal_correction.csv (g_thermal_correction_hartree).  Per state, G_thermal is made
relative to the state's lowest HSE06 value (kcal/mol) and compared: Pearson r, Spearman rho,
MAE and max |error| of the xtb relative G against HSE06.  Acceptance (user rule, 2026-10-21):
max |error| of the order of 0.3-0.5 kcal/mol or smaller.
Conformers whose bhess.out reports a significant imaginary frequency are discarded (both variants).

Usage:
  python analyze_xtb_gtherm_pilot.py --root ~/work_turbomole/artepillin_C_studies --out paper_2/Figures/data/xtb_gtherm_pilot.csv
"""
import argparse
import csv
import os
import re

import numpy as np

HARTREE_TO_KCAL = 627.509474
STATES = ("neutral", "monoanion_a", "monoanion_b", "dianion")
G_RRHO = re.compile(r"G\(RRHO\) contrib\.\s+(-?\d+\.\d+)\s+Eh")


def parse_g_rrho(path):
    """G(RRHO) contribution (Eh) from an xtb output, or None."""
    if not os.path.isfile(path):
        return None
    m = G_RRHO.findall(open(path, errors="ignore").read())
    return float(m[-1]) if m else None


def has_imaginary(root, state, job):
    """True if xtb bhess found a significant imaginary frequency (not a stationary point: discarded)."""
    path = os.path.join(root, "xtb_gtherm_pilot", state, job, "bhess.out")
    return os.path.isfile(path) and "significant imaginary frequency" in open(path, errors="ignore").read()


def relative_kcal(values):
    v = np.asarray(values, float)
    return (v - v.min()) * HARTREE_TO_KCAL


def spearman(a, b):
    ra = np.argsort(np.argsort(a))
    rb = np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


def compare(hse, xtb):
    """hse, xtb: absolute G_thermal in Eh for the same conformers. Returns metrics on the relative scale."""
    h, x = relative_kcal(hse), relative_kcal(xtb)
    # reference both to the HSE06 lowest conformer so the offset cancels
    x = (np.asarray(xtb, float) - np.asarray(xtb, float)[int(np.argmin(hse))]) * HARTREE_TO_KCAL
    err = x - h
    out = dict(n=len(h), mae=float(np.abs(err).mean()), max_abs=float(np.abs(err).max()))
    if len(h) > 2:
        out["pearson"] = float(np.corrcoef(h, x)[0, 1])
        out["spearman"] = spearman(h, x)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    rows = []
    for state in STATES:
        with open(os.path.join(root, f"artepillin_c_{state}_d3_hse06", "layer_b", "freeh_thermal_correction.csv")) as fh:
            hse = {r["job_dir"]: float(r["g_thermal_correction_hartree"]) for r in csv.DictReader(fh)}
        for variant in ("bhess", "hess"):
            jobs = [j for j in hse if not has_imaginary(root, state, j) and parse_g_rrho(os.path.join(root, "xtb_gtherm_pilot", state, j, f"{variant}.out")) is not None]
            if len(jobs) < 2:
                print(f"{state} {variant}: {len(jobs)} usable conformers, skipped")
                continue
            xt = [parse_g_rrho(os.path.join(root, "xtb_gtherm_pilot", state, j, f"{variant}.out")) for j in jobs]
            m = compare([hse[j] for j in jobs], xt)
            rows.append(dict(state=state, variant=variant, **{k: round(v, 3) if isinstance(v, float) else v for k, v in m.items()}))
            print(rows[-1])
    if rows:
        keys = sorted({k for r in rows for k in r}, key=lambda k: ["state", "variant", "n", "mae", "max_abs", "pearson", "spearman"].index(k))
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
