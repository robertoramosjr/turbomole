#!/usr/bin/env python3
"""
analyze_alpb_rerank.py

Compares the VACUUM GFN2 CREST ranking with the GFN2/ALPB(water)
re-optimized ranking of one protonation state
(protonation_and_conformer_investigation.md Sec. 9.38; jobs from
job_xtb_alpb_rerank.sh). Reports, per state:
  * Boltzmann population (298.15 K, CREST degeneracy from
    layer_a/population_table.csv) in vacuum vs in water;
  * which conformers were kept by the original 95% vacuum cut and how
    much of the WATER population they hold;
  * radius of gyration of the population-weighted ensemble in both
    (folding indicator);
  * the conformers that would enter/leave a 95% water cut.

Usage:
    python analyze_alpb_rerank.py <state_dir> [--cut 0.95]
"""

import argparse
import csv
import math
import os
import re

import numpy as np

KCAL = 627.509474
RT = 0.0019872 * 298.15


def read_xyz(p):
    with open(p) as f:
        n = int(f.readline())
        f.readline()
        a = [l.split() for l in f.readlines()[:n]]
    return [x[0] for x in a], np.array([[float(v) for v in x[1:4]] for x in a])


def rg(xyz):
    _, c = read_xyz(xyz)
    return math.sqrt(np.mean(np.sum((c - c.mean(0)) ** 2, axis=1)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("state_dir")
    ap.add_argument("--cut", type=float, default=0.95)
    args = ap.parse_args()
    rows = list(csv.DictReader(open(os.path.join(args.state_dir, "layer_a", "population_table.csv"))))
    rec = []
    for r in rows:
        i = int(r["crest_index"])
        d = os.path.join(args.state_dir, "alpb_rerank", "conf_%03d" % i)
        out = os.path.join(d, "xtb.out")
        e = None
        if os.path.isfile(out):
            m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", open(out).read())
            if m and os.path.isfile(os.path.join(d, "xtbopt.xyz")):
                e = float(m[-1])
        rec.append({"i": i, "deg": int(r["degeneracy"]), "kept": r["kept"] == "1",
                    "ev": float(r["relative_energy_kcal"]), "ew": e, "d": d})
    missing = [r["i"] for r in rec if r["ew"] is None]
    if missing:
        print("WARNING: no ALPB result for conformers", missing)
    rec = [r for r in rec if r["ew"] is not None]
    e0 = min(r["ew"] for r in rec)
    for r in rec:
        r["ewr"] = (r["ew"] - e0) * KCAL
    def pops(key):
        w = np.array([r["deg"] * math.exp(-r[key] / RT) for r in rec])
        return w / w.sum()
    pv, pw = pops("ev"), pops("ewr")
    for r, a, b in zip(rec, pv, pw):
        r["pv"], r["pw"] = a, b
    print(f"{'idx':>4s} {'deg':>4s} {'kept':>4s} {'E_vac':>7s} {'E_alpb':>7s} {'p_vac':>7s} {'p_alpb':>7s} {'Rg_vac':>7s} {'Rg_alpb':>8s}")
    for r in sorted(rec, key=lambda r: -r["pw"])[:25]:
        rv = rg(os.path.join(r["d"], "start.xyz"))
        ra = rg(os.path.join(r["d"], "xtbopt.xyz"))
        r["rgv"], r["rga"] = rv, ra
        print(f"{r['i']:4d} {r['deg']:4d} {int(r['kept']):4d} {r['ev']:7.2f} {r['ewr']:7.2f} {r['pv']:7.3f} {r['pw']:7.3f} {rv:7.3f} {ra:8.3f}")
    for r in rec:
        r.setdefault("rgv", rg(os.path.join(r["d"], "start.xyz")))
        r.setdefault("rga", rg(os.path.join(r["d"], "xtbopt.xyz")))
    kept_w = sum(r["pw"] for r in rec if r["kept"])
    print(f"\nwater population held by the conformers kept by the vacuum {args.cut:.0%} cut: {kept_w:.3f}")
    order = sorted(rec, key=lambda r: -r["pw"]); c = 0; n95 = 0
    for r in order:
        c += r["pw"]; n95 += 1
        if c >= args.cut: break
    new = [r["i"] for r in order[:n95] if not r["kept"]]
    print(f"conformers needed for {args.cut:.0%} of the WATER population: {n95}; not in the original kept set: {new}")
    print("Rg (A), pop-weighted:  vacuum-ensemble/start geoms %.3f | water ensemble/ALPB geoms %.3f | water-pop on start geoms %.3f" % (
        sum(r["pv"] * r["rgv"] for r in rec), sum(r["pw"] * r["rga"] for r in rec), sum(r["pw"] * r["rgv"] for r in rec)))
    sh = np.array([abs(a - b) for a, b in zip(pv, pw)]).sum() / 2
    print(f"total-variation distance between vacuum and water populations: {sh:.3f}")


if __name__ == "__main__":
    main()
