#!/usr/bin/env python3
"""
build_ferulate_reference.py

Builds trans-ferulate (4-hydroxy-3-methoxycinnamate) monoanion (carboxylate + phenol OH) and dianion (carboxylate +
phenolate) starting geometries for the pKa reference pipeline, from the hand-built planar 4-coumarate xyz files
(the 3-position H is replaced by OCH3). Candidate rotamers (OMe methyl pointing toward / away from the phenol O;
phenol H toward / away from OMe, monoanion only) are all relaxed with GFN2-xTB and the lowest-energy one is written.

Usage:
    python build_ferulate_reference.py --coum-mono coum_mono.xyz --coum-di coum_di.xyz \
        --xtb ~/software/TURBOMOLE/bin/em64t-unknown-linux-gnu/xtb --workdir ferulate_build --outdir .
Outputs: ferul_mono.xyz, ferul_di.xyz (+ a text report of the candidate energies).
"""
import argparse
import os
import re
import subprocess

import numpy as np


def read_xyz(p):
    L = open(p).read().split("\n")
    n = int(L[0])
    return L[1], [(l.split()[0], np.array([float(x) for x in l.split()[1:4]])) for l in L[2:2 + n]]


def write_xyz(p, title, atoms):
    with open(p, "w") as f:
        f.write(f"{len(atoms)}\n{title}\n")
        for el, c in atoms:
            f.write(f"{el:2s} {c[0]:12.6f} {c[1]:12.6f} {c[2]:12.6f}\n")


def rot(v, ang):
    c, s = np.cos(ang), np.sin(ang)
    return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1], 0.0])


def add_ome(atoms, ome_sign, oh_toward):
    """Replace H (index 6, ortho to the phenol O at index 10) by OCH3; set phenol H orientation if present."""
    at = [(e, c.copy()) for e, c in atoms]
    c_ring, h = at[1][1], at[6][1]
    u = (h - c_ring) / np.linalg.norm(h - c_ring)
    o = c_ring + 1.36 * u
    back = (c_ring - o) / np.linalg.norm(c_ring - o)
    d = rot(back, ome_sign * np.radians(117.0))
    cm = o + 1.43 * d
    a = d
    uu = np.array([-a[1], a[0], 0.0])  # in-plane perpendicular
    vv = np.array([0.0, 0.0, 1.0])
    hs = [cm + 1.09 * (np.cos(np.radians(70.5)) * a + np.sin(np.radians(70.5)) * (np.cos(p) * uu + np.sin(p) * vv))
          for p in np.radians([0, 120, 240])]
    del at[6]
    at.insert(6, ("O", o))
    at.append(("C", cm))
    at += [("H", x) for x in hs]
    # phenol H (last original atom, if it is an H bonded to O at index 10/ now index 10)
    for i, (e, c) in enumerate(at):
        if e == "H" and np.linalg.norm(c - at[10][1]) < 1.1:
            op = at[10][1]
            r0 = at[0][1]
            base = (op - r0) / np.linalg.norm(op - r0)
            toward = rot(base, np.radians(70.0) * (1 if oh_toward else -1)) if False else None
            sgn = 1 if oh_toward else -1  # +y side = side of the OMe at C(index 1)
            at[i] = ("H", op + 0.96 * rot(base, sgn * np.radians(70.0)))
    return at


def xtb_energy(xtb, xyz, chrg, wd):
    r = subprocess.run([xtb, xyz, "--gfn", "2", "--opt", "tight", "--chrg", str(chrg), "--alpb", "water"],
                       cwd=wd, capture_output=True, text=True).stdout
    m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", r)
    return float(m[-1]) if m else None


def main():
    ap = argparse.ArgumentParser(description="Build ferulate reference geometries.")
    ap.add_argument("--coum-mono", required=True)
    ap.add_argument("--coum-di", required=True)
    ap.add_argument("--xtb", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    rep = []
    for tag, src, chrg in (("mono", a.coum_mono, -1), ("di", a.coum_di, -2)):
        _, atoms = read_xyz(src)
        best = None
        for sign in (1, -1):
            for toward in ((True, False) if tag == "mono" else (True,)):
                cand = add_ome(atoms, sign, toward)
                name = f"{tag}_s{sign}_t{int(toward)}"
                d = os.path.join(a.workdir, name)
                os.makedirs(d, exist_ok=True)
                write_xyz(os.path.join(d, "in.xyz"), name, cand)
                e = xtb_energy(a.xtb, "in.xyz", chrg, d)
                rep.append(f"{name}  E(GFN2/ALPB)={e}")
                if e is not None and (best is None or e < best[0]):
                    best = (e, os.path.join(d, "xtbopt.xyz"), name)
        _, opt = read_xyz(best[1])
        write_xyz(os.path.join(a.outdir, f"ferul_{tag}.xyz"),
                  f"trans-ferulate {'monoanion' if tag == 'mono' else 'dianion'}, GFN2/ALPB best rotamer {best[2]}", opt)
        rep.append(f"-> {tag}: best {best[2]}")
    open(os.path.join(a.outdir, "ferulate_build_report.txt"), "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
