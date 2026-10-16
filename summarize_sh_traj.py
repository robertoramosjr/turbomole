#!/usr/bin/env python3
"""
summarize_sh_traj.py

Per-trajectory status of a native-TURBOMOLE surface-hopping run (phase 1 = S1 SH, phase 2 = S0 continuation).
status: hop_ok | no_hop | died_S1 | died_S0.  product (cis|trans) only when status == hop_ok.
Reuses read()/dihedral() from classify_tm_md.py.

Usage: python3 summarize_sh_traj.py <traj_dir> [...] [--csv out.csv] [--s0-dir DIR]
"""
import argparse
import csv
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classify_tm_md import read  # noqa: E402

DT_FS = 0.5
MIN_S0_FS = 20.0
LAST_FS = 10.0
ATOMS = [0, 22, 23, 6]
NAT = 24
COLS = ["traj", "status", "t_hop_fs", "cnnc_hop", "n_hops", "t_end_s1_fs", "t_end_s0_fs",
        "cnnc_mean_last10", "product"]


def _hop_steps(path):
    if not os.path.exists(path):
        return []
    txt = open(path).read()
    steps, cur = [], False
    for line in txt.split("\n"):
        if re.match(r"\s*-\s*event:\s*hop", line):
            cur = True
        elif cur and (m := re.match(r"\s*steps:\s*(\d+)", line)):
            steps.append(int(m.group(1)))
            cur = False
    return steps


def _died_abnormally(jobex_out):
    if not os.path.exists(jobex_out):
        return False
    tail = [l for l in open(jobex_out).read().split("\n") if l.strip()][-8:]
    return any("ended abnormally" in l for l in tail)


def summarize(traj_dir, s0_dir=None):
    traj_dir = traj_dir.rstrip("/")
    s0_dir = s0_dir or os.path.join(traj_dir, "s0")
    r = dict.fromkeys(COLS)
    r["traj"] = os.path.basename(traj_dir)
    f1 = read(os.path.join(traj_dir, "mdlog.1"), NAT, ATOMS) if os.path.exists(
        os.path.join(traj_dir, "mdlog.1")) else []
    r["t_end_s1_fs"] = f1[-1][0] if f1 else None
    hops = _hop_steps(os.path.join(traj_dir, "event.yaml"))
    r["n_hops"] = len(hops)
    if not hops:
        r["status"] = "died_S1" if _died_abnormally(os.path.join(traj_dir, "jobex.out")) else "no_hop"
        return r
    r["t_hop_fs"] = hops[0] * DT_FS
    if f1:
        r["cnnc_hop"] = min(f1, key=lambda x: abs(x[0] - r["t_hop_fs"]))[1]
    p0 = os.path.join(s0_dir, "mdlog.1")
    f0 = read(p0, NAT, ATOMS) if os.path.exists(p0) else []
    r["t_end_s0_fs"] = f0[-1][0] if f0 else None
    if not f0 or f0[-1][0] - r["t_hop_fs"] < MIN_S0_FS:
        r["status"] = "died_S0"
        return r
    t_end = f0[-1][0]
    m = float(np.mean([abs(x[1]) for x in f0 if x[0] >= t_end - LAST_FS]))
    r["cnnc_mean_last10"] = m
    r["status"] = "hop_ok"
    r["product"] = "trans" if m >= 90 else "cis"
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--csv")
    ap.add_argument("--s0-dir", help="only meaningful with a single traj dir")
    a = ap.parse_args()
    rows = [summarize(d, a.s0_dir) for d in a.dirs]
    for r in rows:
        print(" ".join(f"{k}={'-' if r[k] is None else (round(r[k], 1) if isinstance(r[k], float) else r[k])}"
                       for k in COLS))
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS)
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
