#!/usr/bin/env python3
"""Absolute epsilon (calc vs Camuri 2018 exp) for Paper 2, plan E1.

eps_max = 1.3063e8 * max(sum_i pop_i f_i exp(-(E-E_i)^2/2s^2)) / (s * 8065.54)   [M-1 cm-1, s in eV]
Run with: LC_ALL=C ~/.conda/envs/veusz_env/bin/python scripts/epsilon_abs_table.py
With --subdir eps_cosmo_tddft the same quantities are computed from the COSMO-TDDFT conformers (anion_band_cosmo.cosmo_jobs),
alongside the vacuum HSE06 values of the SAME conformers/weights (like-for-like), and written to a separate CSV.
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_paper2_figures import HC, sticks, state_spectrum  # noqa: E402
from anion_band_cosmo import cosmo_jobs  # noqa: E402

K_EPS = 1.3063e8   # M-1 cm-1 * cm-1 per unit f (Gaussian std sigma)
EV_CM = 8065.54    # cm-1 per eV
SQ2PI = np.sqrt(2 * np.pi)
EXP = {"N": ("s00", 0.61e4), "A": ("s23", 2.13e4)}


def eps_from_profile(profile, grid_ev, sigma_ev, window_ev):
    m = (grid_ev >= window_ev[0]) & (grid_ev <= window_ev[1])
    return K_EPS * profile[m].max() / (sigma_ev * EV_CM)


def subset_spectrum(paths, weights, grid, sigma):
    """Population-weighted Gaussian sum over the given escf outputs; weights renormalised over the subset."""
    tot = sum(weights)
    prof = np.zeros_like(grid)
    for p, w in zip(paths, weights):
        for e, f in sticks(p):
            prof += w / tot * f * np.exp(-((grid - e) ** 2) / (2 * sigma ** 2))
    return prof


def fwhm_ev(wl_nm, y):
    """FWHM in energy of the peak at max(y), contiguous half-height crossing; (nan, note) if a side is cut."""
    ok = ~np.isnan(y)
    wl, y = wl_nm[ok], y[ok]
    e = HC / wl
    o = np.argsort(e)
    e, y = e[o], y[o]
    i = int(np.argmax(y))
    h = y[i] / 2
    lo = i
    while lo > 0 and y[lo] > h:
        lo -= 1
    hi = i
    while hi < len(y) - 1 and y[hi] > h:
        hi += 1
    cut = y[lo] > h or y[hi] > h
    el = np.interp(h, [y[lo], y[lo + 1]], [e[lo], e[lo + 1]]) if y[lo] <= h else np.nan
    eh = np.interp(h, [y[hi], y[hi - 1]], [e[hi], e[hi - 1]]) if y[hi] <= h else np.nan
    return eh - el, ("half-height not reached on one side" if cut else "")


def load_exp(path):
    rows = [ln.split() for ln in open(path) if not ln.startswith("#")]
    head, data = rows[0], rows[1:]
    arr = np.array([[float(v) for v in r] for r in data])  # 'nan' parses
    return {h: arr[:, k] for k, h in enumerate(head)}


def write_cosmo(a, grid, win, w2, sig_exp):
    """Rows medium in {vacuum, cosmo} over the COSMO-TDDFT conformers; ratio_to_N_same_medium gives the A/N intensity ratio."""
    rows, eps = [], {}
    for s, (col, eexp) in EXP.items():
        jobs, pop = cosmo_jobs(a.root, s)
        w = [pop[j] for j, _, _ in jobs]
        for label, sg in (("paper", a.sigma), ("exp", sig_exp[s])):
            for medium, k in (("vacuum", 1), ("cosmo", 2)):
                prof = subset_spectrum([t[k] for t in jobs], w, grid, sg)
                e = eps_from_profile(prof, grid, sg, w2)
                m = (grid >= win[0]) & (grid <= win[1])
                area = np.trapezoid(prof[m], grid[m]) / (sg * SQ2PI)
                eps[(s, label, medium)] = e
                rows.append(dict(state=s, medium=medium, sigma_eV=f"{sg:.4f}", sigma_basis=label, eps_max_calc=f"{e:.0f}",
                                 eps_exp=f"{eexp:.0f}", ratio_calc_exp=f"{e / eexp:.2f}", n_conf=len(jobs),
                                 coverage_pop=f"{sum(w) / sum(pop.values()):.3f}", sum_pop_f_250_400nm=f"{area:.3f}"))
    for r in rows:
        r["ratio_A_over_N_same_medium"] = f"{eps[('A', r['sigma_basis'], r['medium'])] / eps[('N', r['sigma_basis'], r['medium'])]:.2f}"
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)
    print(f"exp A/N = {2.13 / 0.61:.2f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default="/home/rramos/work_turbomole/artepillin_C_studies")
    ap.add_argument("--exp", default="Figures/data/exp_camuri2018_fig2_spectra_vs_pH.dat")
    ap.add_argument("--out", default=None, help="default depends on --subdir")
    ap.add_argument("--subdir", default=None, help="eps_cosmo_tddft: use the COSMO-TDDFT conformers (vacuum vs COSMO, same jobs)")
    ap.add_argument("--sigma", type=float, default=0.4)
    a = ap.parse_args()
    if a.out is None:
        a.out = "Figures/data/epsilon_calc_vs_exp_cosmo.csv" if a.subdir else "Figures/data/epsilon_calc_vs_exp.csv"

    exp = load_exp(a.exp)
    wl = exp["wavelength_nm"]
    grid = np.linspace(2.0, 7.0, 5001)
    win = (HC / 400, HC / 250)
    w2 = (3.0, 5.5)
    rows = []
    sig_exp = {}
    for s, (col, eexp) in EXP.items():
        fw, note = fwhm_ev(wl, exp[col])
        sig_exp[s] = fw / 2.3548
        print(f"{s}: FWHM_exp = {fw:.3f} eV sigma_exp = {sig_exp[s]:.3f} eV {note}")
    if a.subdir:
        return write_cosmo(a, grid, win, w2, sig_exp)
    area = {}
    for s, (col, eexp) in EXP.items():
        for label, sg in (("paper", a.sigma), ("exp", sig_exp[s])):
            prof, nl, nt, cov = state_spectrum(a.root, s, grid, sg)
            e = eps_from_profile(prof, grid, sg, w2)
            m = (grid >= win[0]) & (grid <= win[1])
            area[(s, label)] = np.trapezoid(prof[m], grid[m]) / (sg * SQ2PI)  # sum pop*f in 250-400 nm
            rows.append(dict(state=s, sigma_eV=f"{sg:.4f}", sigma_basis=label, eps_max_calc=f"{e:.0f}",
                             eps_exp=f"{eexp:.0f}", ratio_calc_exp=f"{e / eexp:.2f}", coverage_pop=f"{cov:.3f}",
                             n_labelled=nl, n_total=nt, fwhm_exp_eV=f"{sig_exp[s] * 2.3548:.3f}",
                             sum_pop_f_250_400nm=f"{area[(s, label)]:.3f}"))
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)
    print(f"A/N: eps(sigma=0.4) = {float(rows[2]['eps_max_calc']) / float(rows[0]['eps_max_calc']):.2f}; "
          f"sum f = {area[('A', 'paper')] / area[('N', 'paper')]:.2f}; exp = {2.13 / 0.61:.2f}")


if __name__ == "__main__":
    main()
