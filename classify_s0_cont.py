#!/usr/bin/env python3
"""Classifica o produto (cis/trans) das continuações em S0 do azobenzeno pelo CNNC final.
Uso: classify_s0_cont.py <.../5-DYN/S0_cont/TRAJECTORIES> [--last-fs 20]
CNNC = C1-N23-N24-C7 (índices 0,22,23,6 no XYZ). |φ|<90 -> cis, >90 -> trans.
Produto = média de |φ| nos últimos --last-fs fs (evita decidir por um frame só)."""
import sys, glob, os, re, argparse
import numpy as np

def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))

def read_xyz(path):
    frames, lines = [], open(path).read().splitlines()
    i = 0
    while i < len(lines):
        n = int(lines[i].split()[0]); t = float(lines[i + 1].split()[0])
        xyz = np.array([[float(x) for x in l.split()[1:4]] for l in lines[i + 2:i + 2 + n]])
        frames.append((t, xyz)); i += n + 2
    return frames

ap = argparse.ArgumentParser(); ap.add_argument("tdir"); ap.add_argument("--last-fs", type=float, default=20)
a = ap.parse_args()
rows = []
for d in sorted(glob.glob(os.path.join(a.tdir, "TRAJ*")), key=lambda s: int(re.findall(r"\d+", s)[-1])):
    n = int(re.findall(r"\d+", d)[-1]); f = os.path.join(d, "RESULTS", "geometries.xyz")
    if not os.path.exists(f): rows.append((n, None, None, None, "SEM_RESULTS")); continue
    fr = read_xyz(f); t = np.array([x[0] for x in fr])
    phi = np.array([dihedral(*[x[1][k] for k in (0, 22, 23, 6)]) for x in fr])
    m = t >= t[-1] - a.last_fs; mean_abs = np.abs(phi[m]).mean()
    rows.append((n, t[-1], phi[-1], mean_abs, "trans" if mean_abs > 90 else "cis"))
ok = [r for r in rows if r[1] is not None]
print("traj  t_final(fs)  CNNC_final  <|CNNC|>_last  produto")
for r in rows:
    print(f"{r[0]:4d}" + (f"  {r[1]:9.1f}  {r[2]:9.1f}  {r[3]:11.1f}  {r[4]}" if r[1] is not None else "  " + r[4]))
nt = sum(r[4] == "trans" for r in ok); nc = sum(r[4] == "cis" for r in ok); N = nt + nc
p = nt / N if N else float("nan"); se = (p * (1 - p) / N) ** 0.5 if N else float("nan")
print(f"\nN={N}  trans={nt}  cis={nc}  fração trans={p:.2f} ± {se:.2f} (erro binomial)")
