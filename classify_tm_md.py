#!/usr/bin/env python3
"""
classify_tm_md.py

Reads one or more TURBOMOLE MD logs (mdlog.1) and reports the CNNC dihedral at the end and its mean over the
last N fs, classifying |phi| < 90 deg as cis-like and >= 90 deg as trans-like (azobenzene; CNNC atoms are
0-based indices 0, 22, 23, 6 of the cis geometry order used in this project, as in classify_s0_cont.py).

Usage:
    python classify_tm_md.py <mdlog.1> [<mdlog.1> ...] [--last-fs 10] [--atoms 0 22 23 6] [--natoms 24]
"""
import argparse
import math

import numpy as np

FS = 41.341373


def dihedral(p):
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    return math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def read(path, nat, atoms):
    lines = open(path).read().split("\n")
    out = []
    for i, l in enumerate(lines):
        if l.startswith("t="):
            try:
                c = np.array([[float(x) for x in lines[i + 1 + k].split()[:3]] for k in range(nat)])
                out.append((float(l[2:]) / FS, dihedral(c[atoms])))
            except (ValueError, IndexError):
                pass
    return out


def main():
    ap = argparse.ArgumentParser(description="CNNC-based product classification of TM MD logs.")
    ap.add_argument("logs", nargs="+")
    ap.add_argument("--last-fs", type=float, default=10.0)
    ap.add_argument("--atoms", type=int, nargs=4, default=[0, 22, 23, 6])
    ap.add_argument("--natoms", type=int, default=24)
    a = ap.parse_args()
    for p in a.logs:
        d = read(p, a.natoms, a.atoms)
        if not d:
            print(p, "no frames")
            continue
        t_end = d[-1][0]
        tail = [abs(x[1]) for x in d if x[0] >= t_end - a.last_fs]
        m = float(np.mean(tail))
        print(f"{p}: t_end={t_end:.1f} fs  |CNNC|_end={abs(d[-1][1]):.0f}  mean|CNNC|(last {a.last_fs:.0f} fs)={m:.0f}  -> "
              f"{'trans-like' if m >= 90 else 'cis-like'}")


if __name__ == "__main__":
    main()
