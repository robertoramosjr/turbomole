#!/usr/bin/env python3
"""
plot_vacuum_vs_water_ensembles.py

Paper figure/table (Sec. 9.39): vacuum-GFN2 vs GFN2/ALPB(water)
conformer populations for the anionic states (dianion, monoanion B,
monoanion A), from job_xtb_alpb_rerank.sh outputs. Panels: cumulative
population vs conformer rank in both media (with the size of the
original 95%-vacuum set marked), and per-conformer population vac vs
water (log-log, kept vs dropped conformers).

Usage:
    python plot_vacuum_vs_water_ensembles.py --project-root . \
        --output-prefix vacuum_vs_water_ensembles
"""

import argparse
import csv
import math
import os
import re

import numpy as np
import matplotlib.pyplot as plt

KCAL = 627.509474
RT = 0.0019872 * 298.15
STATES = [("dianion", "Dianion"), ("monoanion_b", "Monoanion B (phenolate)"), ("monoanion_a", "Monoanion A (carboxylate)")]


def load(state_dir):
    rows = list(csv.DictReader(open(os.path.join(state_dir, "layer_a", "population_table.csv"))))
    rec = []
    for r in rows:
        i = int(r["crest_index"])
        out = os.path.join(state_dir, "alpb_rerank", "conf_%03d" % i, "xtb.out")
        m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", open(out).read()) if os.path.isfile(out) else []
        if m:
            rec.append([i, int(r["degeneracy"]), r["kept"] == "1", float(r["relative_energy_kcal"]), float(m[-1])])
    e0 = min(r[4] for r in rec)
    a = np.array([[r[1], r[3], (r[4] - e0) * KCAL] for r in rec])
    kept = np.array([r[2] for r in rec])
    pv = a[:, 0] * np.exp(-a[:, 1] / RT); pv /= pv.sum()
    pw = a[:, 0] * np.exp(-a[:, 2] / RT); pw /= pw.sum()
    return pv, pw, kept


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--output-prefix", required=True)
    args = ap.parse_args()
    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    lines = ["state,n_crest,n_kept_vac95,water_pop_in_kept,n_for_95_vac,n_for_95_water,tv_distance"]
    for c, (s, title) in enumerate(STATES):
        pv, pw, kept = load(os.path.join(args.project_root, f"artepillin_c_{s}_d3_hse06"))
        cv = np.cumsum(np.sort(pv)[::-1]); cw = np.cumsum(np.sort(pw)[::-1])
        n95v = int(np.searchsorted(cv, 0.95) + 1); n95w = int(np.searchsorted(cw, 0.95) + 1)
        x = np.arange(1, len(pv) + 1)
        ax[0, c].plot(x, cv, label="vacuum GFN2"); ax[0, c].plot(x, cw, label="GFN2/ALPB(water)")
        ax[0, c].axhline(0.95, color="grey", ls=":"); ax[0, c].axvline(kept.sum(), color="k", ls="--", lw=0.8,
                          label=f"original set ({int(kept.sum())})")
        ax[0, c].set_xscale("log"); ax[0, c].set_xlabel("conformer rank"); ax[0, c].set_ylabel("cumulative population")
        ax[0, c].set_title(f"{title}\n95% needs {n95v} (vacuum) / {n95w} (water)", fontsize=10); ax[0, c].legend(fontsize=8, loc="lower right")
        ax[1, c].scatter(pv[kept], pw[kept], s=14, label="in original set", alpha=.8)
        ax[1, c].scatter(pv[~kept], pw[~kept], s=14, marker="x", label="dropped by vacuum cut", alpha=.8)
        lim = [max(min(pv.min(), pw.min()), 1e-6), 1]; ax[1, c].plot(lim, lim, "k-", lw=.5)
        ax[1, c].set_xscale("log"); ax[1, c].set_yscale("log"); ax[1, c].set_xlim(lim); ax[1, c].set_ylim(lim)
        ax[1, c].set_xlabel("population, vacuum"); ax[1, c].set_ylabel("population, water"); ax[1, c].legend(fontsize=8)
        tv = 0.5 * np.abs(pv - pw).sum()
        lines.append(f"{s},{len(pv)},{int(kept.sum())},{pw[kept].sum():.3f},{n95v},{n95w},{tv:.3f}")
    fig.tight_layout(); fig.savefig(args.output_prefix + ".png", dpi=200)
    with open(args.output_prefix + ".dat", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
