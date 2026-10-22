#!/usr/bin/env python3
"""Prepare COSMO-TDDFT (escf) runs for the aqueous conformers selected by the free-energy weights.

Per state the conformers with the largest weight w_G (water_gens_conformers.csv from collect_water_gens.py)
are taken until their cumulative weight reaches --cover.  Each W4 directory (HSE06-D3(BJ)/def2-TZVP + COSMO
ridft, converged) is copied to <root>/water_tddft/<state>/<cNNNN>/ (originals untouched) and the escf
keywords of the production protocol are appended to `control`: $scfinstab rpas, $soes a <roots>, $rpacor 500
(fast response through the COSMO refractive index, as in the eps_cosmo_tddft runs).
The joblist lists the top-weighted conformer of each state first (pilot: array 1-4), then the rest.

Usage:
  python prepare_water_escf.py --root ~/work_turbomole/artepillin_C_studies --cover 0.95 [--roots 60] [--dry-run]
"""
import argparse
import csv
import os
import re
import shutil

STATES = {"N": "neutral", "A": "monoanion_a", "B": "monoanion_b", "D": "dianion"}
KEEP = ("control", "coord", "basis", "auxbasis", "mos", "out.ccf", "cosmo_transfer.tmp")


def select(rows, cover):
    """rows: [(job, weight)] -> jobs by decreasing weight until cumulative >= cover."""
    tot = sum(w for _, w in rows)
    acc, sel = 0.0, []
    for job, w in sorted(rows, key=lambda r: -r[1]):
        sel.append(job)
        acc += w / tot
        if acc >= cover:
            break
    return sel


def prepare(src, dst, roots):
    os.makedirs(dst, exist_ok=True)
    for f in KEEP:
        if os.path.isfile(os.path.join(src, f)):
            shutil.copy(os.path.join(src, f), dst)
    ctl = open(os.path.join(dst, "control")).read()
    if re.search(r"^\$soes", ctl, re.M) or "$end" not in ctl:
        raise SystemExit(f"{dst}: control already has $soes or lacks $end")
    block = f"$scfinstab rpas\n$soes\n a {roots:>11d}\n$rpacor 500\n"
    open(os.path.join(dst, "control"), "w").write(ctl.rsplit("$end", 1)[0] + block + "$end\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--cover", type=float, default=0.95)
    ap.add_argument("--roots", type=int, default=60)
    ap.add_argument("--weights", default="paper_2/Figures/data/water_gens_conformers.csv")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    data = list(csv.DictReader(open(os.path.join(root, a.weights))))
    first, rest = [], []
    for label, state in STATES.items():
        rows = [(r["job"], float(r["w_G"])) for r in data if r["state"] == label and r["w_G"]]
        sel = select(rows, a.cover)
        print(f"{label}: {len(sel)} of {len(rows)} conformers (cover {a.cover})")
        for k, job in enumerate(sel):
            src = os.path.join(root, "w4_sp", state, job)
            dst = os.path.join(root, "water_tddft", state, job)
            if not a.dry_run:
                prepare(src, dst, a.roots)
            (first if k == 0 else rest).append(dst)
    allp = first + rest
    if not a.dry_run:
        with open(os.path.join(root, "water_tddft", "joblist.txt"), "w") as f:
            f.write("\n".join(allp) + "\n")
    print(len(allp), "jobs; first", len(first), "= pilot (one per state)")


if __name__ == "__main__":
    main()
