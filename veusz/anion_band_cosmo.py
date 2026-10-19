"""H4: lambda_max of N and A with COSMO(water) vs vacuum, same HSE06/def2-TZVP conformers and weights.

For each state, the conformers in eps_cosmo_tddft/<N|A>/job_NNNN are compared with their vacuum HSE06 escf
(layer_b); Boltzmann weights (r2scan3c_population) renormalised over the conformers computed in both.
Writes Figures/data/anion_band_cosmo.csv (state, medium, n_conf, pop_cov, lambda_max_nm, shift vs vacuum).
Window defaults to 3.0-4.96 eV (250-400 nm, the experimental band; higher roots would otherwise win at the window edge).
Usage: LC_ALL=C ~/.conda/envs/veusz_env/bin/python anion_band_cosmo.py [--sigma 0.4] [--emin 3.0 --emax 4.96]
"""
import argparse, csv, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from anion_band_coverage import HC, SDIR, read_pop  # noqa: E402
from anion_band_states import escf_path  # noqa: E402
from plot_paper2_figures import sticks  # noqa: E402


def lambda_max(paths, weights, sigma, emin, emax):
    """Weighted Gaussian-sum band maximum (nm) inside [emin, emax] eV; roots outside the window are ignored."""
    grid = np.linspace(emin, emax, 4000)
    tot = sum(weights)
    prof = np.zeros_like(grid)
    for p, w in zip(paths, weights):
        for e, f in sticks(p):
            prof += w / tot * f * np.exp(-((grid - e) ** 2) / (2 * sigma ** 2))
    return float(HC / grid[np.argmax(prof)])


def cosmo_jobs(root, s):
    """Conformers of state s with a finished COSMO escf and a vacuum HSE06 escf: ([(job, vacuum_escf, cosmo_escf)], pop)."""
    d = os.path.join(root, f"artepillin_c_{SDIR[s]}_d3_hse06")
    pop = read_pop(d)
    jobs = []
    for j in sorted(os.listdir(os.path.join(root, "eps_cosmo_tddft", s))):
        c = os.path.join(root, "eps_cosmo_tddft", s, j, "escf.out")
        v = escf_path(d, j)
        if v and os.path.isfile(c) and "escf : all done" in open(c).read() and j in pop:
            jobs.append((j, v, c))
    return jobs, pop


def run(root, sigma, emin, emax):
    rows = []
    for s in SDIR:
        jobs, pop = cosmo_jobs(root, s)
        w = [pop[j] for j, _, _ in jobs]
        res = {}
        for medium, k in (("vacuum", 1), ("cosmo", 2)):
            res[medium] = lambda_max([t[k] for t in jobs], w, sigma, emin, emax)
        for medium in ("vacuum", "cosmo"):
            rows.append(dict(state=s, medium=medium, n_conf=len(jobs), pop_cov=round(sum(w) / sum(pop.values()), 4),
                             lambda_max_nm=round(res[medium], 1), shift_vs_vacuum_nm=round(res[medium] - res["vacuum"], 1)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default="/home/rramos/work_turbomole/artepillin_C_studies")
    ap.add_argument("--out", default="Figures/data/anion_band_cosmo.csv")
    ap.add_argument("--sigma", type=float, default=0.4)
    ap.add_argument("--emin", type=float, default=3.0)
    ap.add_argument("--emax", type=float, default=4.96)
    a = ap.parse_args()
    rows = run(a.root, a.sigma, a.emin, a.emax)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
