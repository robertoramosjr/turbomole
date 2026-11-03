#!/usr/bin/env python3
"""Map the CREST structures kept by the vacuum 95 % cut onto the aqueous (ALPB) CREST search.

Both searches start from the same input, so atom order is shared.  For each aqueous structure j the script
finds the vacuum CREST structure k with the smallest heavy-atom RMSD (Kabsch), minimised over the heavy-atom
automorphism images (ring flip, gem-dimethyl swaps, carboxylate O swap).  The aqueous structure is 'covered'
if its best match is a vacuum structure kept in layer A (layer_a_joblist.txt, job_NNNN = CREST index) and
RMSD <= --tol.  The retention is the aqueous r2SCAN-3c/COSMO population (w_r2, water_ensemble_ranking.csv,
degeneracy weighted) on covered structures.

Usage:
  python map_vacuum_to_water.py --root ~/work_turbomole/artepillin_C_studies --tol 0.5 \
      --out todo_v12/analysis/vacuum_to_water_map.csv --summary todo_v12/analysis/vacuum_retention_r2cosmo.csv
"""
import argparse
import csv
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from degeneracy_vs_permutations import STATES, count_automorphisms, graph, read_coord  # noqa: E402


def read_multixyz(path):
    lines = open(path).read().splitlines()
    out, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            break
        n = int(lines[i])
        el = [l.split()[0].lower() for l in lines[i + 2:i + 2 + n]]
        xyz = np.array([[float(v) for v in l.split()[1:4]] for l in lines[i + 2:i + 2 + n]])
        out.append((el, xyz))
        i += n + 2
    return out


def kabsch_rmsd(a, b):
    a = a - a.mean(0)
    b = b - b.mean(0)
    u, s, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    s[-1] *= d
    return np.sqrt(max((np.sum(a * a) + np.sum(b * b) - 2 * s.sum()) / len(a), 0.0))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--tol", type=float, default=0.5, help="heavy-atom RMSD tolerance (A)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", required=True)
    a = ap.parse_args()
    root = os.path.expanduser(a.root)
    kept = {}
    for line in open(os.path.join(root, "layer_a_joblist.txt")):
        m = re.search(r"artepillin_c_(\w+?)_d3_hse06/layer_a/conformer_jobs/job_(\d+)", line)
        if m:
            kept.setdefault(m.group(1), set()).add(int(m.group(2)))
    vac_dir = {"N": "neutral", "A": "monoanion_a", "B": "monoanion_b", "D": "dianion"}
    rank = list(csv.DictReader(open(os.path.join(root, "paper_2/Figures/data/water_ensemble_ranking.csv"))))
    rows, summ = [], []
    for st, d in STATES.items():
        el0, xyz0 = read_coord(os.path.join(root, "water_ensemble", d, "coord"))
        nb = graph(el0, xyz0)
        imgs = count_automorphisms(el0, nb, only_heavy="images")
        heavy = [i for i, e in enumerate(el0) if e != "h"]
        pos = {h: p for p, h in enumerate(heavy)}
        perms = [np.array([pos[x] for x in t]) for t in imgs]
        vac = read_multixyz(os.path.join(root, f"artepillin_c_{d}_d3_hse06", "crest_conformers.xyz"))
        wat = read_multixyz(os.path.join(root, "water_ensemble", d, "crest_conformers.xyz"))
        assert vac[0][0] == wat[0][0], "atom order differs between vacuum and water searches"
        vh = [v[1][heavy] for v in vac]
        w_r2 = {int(r["crest_index"]): float(r["w_r2"]) for r in rank if r["state"] == st}
        cov = tot = 0.0
        nmatch = 0
        for j, (_, xw) in enumerate(wat, start=1):
            hw = xw[heavy]
            best = (9e9, 0)
            for k, hv in enumerate(vh, start=1):
                r = min(kabsch_rmsd(hw[p], hv) for p in perms)
                if r < best[0]:
                    best = (r, k)
            ok = best[0] <= a.tol and best[1] in kept.get(d, set())
            rows.append({"state": st, "water_index": j, "vac_index": best[1], "rmsd_heavy": round(best[0], 3),
                         "vac_kept": int(best[1] in kept.get(d, set())), "covered": int(ok), "w_r2": w_r2.get(j, 0.0)})
            tot += w_r2.get(j, 0.0)
            cov += w_r2.get(j, 0.0) * ok
            nmatch += best[0] <= a.tol
        summ.append({"state": st, "tol_A": a.tol, "n_water": len(wat), "n_vac": len(vac), "n_vac_kept": len(kept.get(d, ())),
                     "n_water_matched": nmatch, "retention_w_r2": round(cov / tot, 4)})
        print(summ[-1], flush=True)
    for path, data in ((a.out, rows), (a.summary, summ)):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)


if __name__ == "__main__":
    main()
