#!/usr/bin/env python3
"""Prepare Raman (egrad polly + intense) copies of finished PBE0/def2-SVP aoforce jobs.

For each state, the conformers that carry --cover of the r2SCAN-3c Boltzmann population (and at most
--max-n) are copied from <state-dir>/layer_a/pbe0_jobs/job_NNNN to <out>/<label>/job_NNNN, keeping only what
egrad/intense need (control, coord, basis, auxbasis, mos, hessian, dipgrad, vib_normal_modes, vibspectrum;
the 11 MB `sing_a` of the TD-DFT step is not copied). The copy's control loses the TD-DFT blocks
(`$scfinstab rpas`, `$soes` + its line, `$rpacor`); `$scfinstab polly` is added by job_raman_polly.sh.
A joblist (absolute paths) and a cost summary are written; nothing is submitted.

Usage:
    python prepare_raman_jobs.py --root ~/work_turbomole/artepillin_C_studies \\
        --state N=artepillin_c_neutral_d3_hse06 --state A=artepillin_c_monoanion_a_d3_hse06 \\
        --cover 0.9 --out raman [--max-n 60] [--dry-run]
    python prepare_raman_jobs.py ... --pilot N      # only the top-weighted conformer of state N
"""
import argparse
import csv
import os
import shutil

KEEP = ("control", "coord", "basis", "auxbasis", "mos", "hessian", "dipgrad", "vib_normal_modes", "vibspectrum")
AOFORCE_CORE_H = 16.0   # 8 cores x ~2 h wall for the top neutral conformer (aoforce.out, 2026-08-29): scale reference


def strip_tddft(lines):
    """Drop `$scfinstab ...`, `$soes` (+ its data lines) and `$rpacor` from control lines; keep everything else."""
    out, skip = [], False
    for ln in lines:
        if ln.startswith("$"):
            skip = ln.startswith("$soes")
            if skip or ln.startswith("$scfinstab") or ln.startswith("$rpacor"):
                continue
        elif skip:
            continue
        out.append(ln)
    return out


def select(pop, have, cover, max_n):
    tot = sum(pop[j] for j in have)
    sel, acc = [], 0.0
    for j in sorted(have, key=lambda j: -pop[j]):
        sel.append(j)
        acc += pop[j] / tot
        if acc >= cover or len(sel) == max_n:
            break
    return sel, acc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--state", action="append", required=True, help="LABEL=state directory (relative to --root)")
    ap.add_argument("--cover", type=float, default=0.9)
    ap.add_argument("--max-n", type=int, default=10 ** 6)
    ap.add_argument("--out", default="raman")
    ap.add_argument("--pilot", default=None, help="state label: select only its top-weighted conformer")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    jobs = []
    for spec in a.state:
        label, sd = spec.split("=")
        if a.pilot and label != a.pilot:
            continue
        d = os.path.join(root, sd, "layer_a")
        with open(os.path.join(d, "population_table_dft.csv")) as f:
            pop = {r["job_dir"]: float(r["r2scan3c_population"]) for r in csv.DictReader(f)}
        have = {j for j in pop if os.path.isfile(os.path.join(d, "pbe0_jobs", j, "hessian"))}
        sel, acc = select(pop, have, a.cover, 1 if a.pilot else a.max_n)
        print(f"{label}: {len(sel)} of {len(have)} conformers with hessian, {acc * 100:.1f}% of the population")
        for j in sel:
            dst = os.path.join(root, a.out, label, j)
            jobs.append(dst)
            if a.dry_run:
                continue
            os.makedirs(dst, exist_ok=True)
            src = os.path.join(d, "pbe0_jobs", j)
            for fn in KEEP:
                shutil.copy2(os.path.join(src, fn), os.path.join(dst, fn))
            with open(os.path.join(dst, "control")) as f:
                lines = f.read().split("\n")
            with open(os.path.join(dst, "control"), "w") as f:
                f.write("\n".join(strip_tddft(lines)))
    print(f"total {len(jobs)} jobs; upper-bound cost if egrad ~ aoforce: {len(jobs) * AOFORCE_CORE_H:.0f} core-h")
    if not a.dry_run and jobs:
        lst = os.path.join(root, a.out, "joblist_pilot.txt" if a.pilot else "joblist.txt")
        with open(lst, "w") as f:
            f.write("\n".join(jobs) + "\n")
        print("joblist:", lst)


if __name__ == "__main__":
    main()
