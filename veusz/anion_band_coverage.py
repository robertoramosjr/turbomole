"""H1: does the HSE06-labelled subset of the ensemble bias the band maximum?

For N and A, compute the PBE0/def2-SVP lambda_max (Gaussian sum in energy) of the
full ensemble and of the subset that also has an HSE06 label, same weights as
state_spectrum (renormalised over the subset).
Run with: LC_ALL=C ~/.conda/envs/veusz_env/bin/python anion_band_coverage.py
"""
import argparse, csv, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_paper2_figures import sticks  # noqa: E402

HC = 1239.84198
SDIR = {"N": "neutral", "A": "monoanion_a"}


def read_pop(d):
    with open(os.path.join(d, "layer_a", "population_table_dft.csv")) as f:
        return {r["job_dir"]: float(r["r2scan3c_population"]) for r in csv.DictReader(f)}


def hse_jobs(d):
    found = set()
    for sub in ("labeled_jobs", "confirmation_batch_jobs", "random_control_batch_jobs"):
        base = os.path.join(d, "layer_b", sub)
        if os.path.isdir(base):
            found |= {j for j in os.listdir(base) if os.path.isfile(os.path.join(base, j, "escf.out"))}
    return found


def lambda_max(jobs, pop, d, sigma, grid):
    tot = sum(pop[j] for j in jobs)
    prof = np.zeros_like(grid)
    for j in jobs:
        p = os.path.join(d, "layer_a", "pbe0_jobs", j, "escf.out")
        for e, f in sticks(p):
            prof += pop[j] / tot * f * np.exp(-((grid - e) ** 2) / (2 * sigma ** 2))
    w = (grid >= 3.0) & (grid <= 5.5)
    emax = float(np.mean([max(e for e, _ in sticks(os.path.join(d, 'layer_a', 'pbe0_jobs', j, 'escf.out'))) for j in jobs]))
    return float(HC / grid[w][np.argmax(prof[w])]), tot, emax


def run(root, sigma):
    grid = np.linspace(2.0, 6.0, 4000)
    rows = []
    for s, name in SDIR.items():
        d = os.path.join(root, f"artepillin_c_{name}_d3_hse06")
        pop = read_pop(d)
        full = {j for j in pop if os.path.isfile(os.path.join(d, "layer_a", "pbe0_jobs", j, "escf.out"))}
        sub = full & hse_jobs(d)
        for lab, jobs in (("full", full), ("hse06_subset", sub)):
            lam, tot, emax = lambda_max(jobs, pop, d, sigma, grid)
            rows.append(dict(state=s, level="PBE0/def2-SVP", subset=lab, n_conf=len(jobs),
                             pop_cov=round(tot / sum(pop.values()), 4), lambda_max_nm=round(lam, 1), e_max_root_ev=round(emax, 2)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--root", default=os.path.abspath(os.path.join(here, "..", "..")))
    ap.add_argument("--sigma", type=float, default=0.4)
    ap.add_argument("--out", default=os.path.join(here, "..", "Figures", "data", "anion_band_coverage.csv"))
    a = ap.parse_args()
    rows = run(a.root, a.sigma)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    for r in rows:
        print(r)
    for s in SDIR:
        l = {r["subset"]: r["lambda_max_nm"] for r in rows if r["state"] == s}
        print(f"{s}: dcov = {l['full'] - l['hse06_subset']:+.1f} nm")


if __name__ == "__main__":
    main()
