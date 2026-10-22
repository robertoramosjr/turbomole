#!/usr/bin/env python3
"""Replace the two para H atoms of azobenzene by an anchor group (-SH or -NH2).

Reads a TURBOMOLE `coord` (Bohr) of azobenzene, writes an XYZ (Angstrom) of the 4,4'-disubstituted
molecule with idealised group geometries (starting guess only; jobex relaxes it).

Usage:
  add_azo_anchors.py --coord azo_trans_hse06_m4_d3/coord --group sh  --output trans_sh.xyz
  add_azo_anchors.py --coord azo_cis_hse06_m4_d3/coord   --group nh2 --output cis_nh2.xyz
"""
import argparse
from pathlib import Path

import numpy as np

BOHR = 0.529177210903
# bond lengths (A) and angles (deg): C-S 1.77, S-H 1.34, C-S-H 96; C-N 1.40, N-H 1.01, C-N-H 121 (planar)
GROUPS = {"sh": ("S", 1.77, [("H", 1.34, 96.0, 1)]),
          "nh2": ("N", 1.40, [("H", 1.01, 121.0, 1), ("H", 1.01, 121.0, -1)])}


def read_coord(path):
    el, xyz = [], []
    for line in Path(path).read_text().split("\n")[1:]:
        if line.startswith("$"):
            break
        p = line.split()
        xyz.append([float(v) * BOHR for v in p[:3]])
        el.append(p[3].capitalize())
    return el, np.array(xyz)


def para_sites(el, xyz):
    """Return [(para C index, its H index)] for each phenyl ring (ring split by nearest N)."""
    n_idx = [i for i, e in enumerate(el) if e == "N"]
    c_idx = [i for i, e in enumerate(el) if e == "C"]
    h_idx = [i for i, e in enumerate(el) if e == "H"]
    rings = {n: [] for n in n_idx}
    for c in c_idx:
        rings[min(n_idx, key=lambda n: np.linalg.norm(xyz[c] - xyz[n]))].append(c)
    sites = []
    for n, ring in rings.items():
        ipso = min(ring, key=lambda c: np.linalg.norm(xyz[c] - xyz[n]))
        para = max(ring, key=lambda c: np.linalg.norm(xyz[c] - xyz[ipso]))
        h = min(h_idx, key=lambda i: np.linalg.norm(xyz[i] - xyz[para]))
        sites.append((para, h, ring))
    return sites


def substitute(el, xyz, group):
    atom, r_cx, hs = GROUPS[group]
    drop, new_el, new_xyz = set(), [], []
    for para, h, ring in para_sites(el, xyz):
        d = xyz[h] - xyz[para]
        d /= np.linalg.norm(d)
        pts = xyz[ring]
        normal = np.cross(pts[1] - pts[0], pts[2] - pts[0])
        normal /= np.linalg.norm(normal)
        p = np.cross(normal, d)
        x_pos = xyz[para] + r_cx * d
        drop.add(h)
        new_el.append(atom)
        new_xyz.append(x_pos)
        for sym, r, ang, sign in hs:
            a = np.radians(ang)
            w = np.cos(a) * (-d) + sign * np.sin(a) * p
            new_el.append(sym)
            new_xyz.append(x_pos + r * w)
    keep = [i for i in range(len(el)) if i not in drop]
    return [el[i] for i in keep] + new_el, np.vstack([xyz[keep], new_xyz])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--coord", required=True, type=Path)
    ap.add_argument("--group", required=True, choices=sorted(GROUPS))
    ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    el, xyz = substitute(*read_coord(a.coord), a.group)
    with open(a.output, "w") as f:
        f.write(f"{len(el)}\n{a.coord.parent.name} + para {a.group.upper()} (starting guess)\n")
        for s, p in zip(el, xyz):
            f.write(f"{s:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")
    heavy = {s: el.count(s) for s in sorted(set(el))}
    print(f"{a.output}: {len(el)} atoms {heavy}")


if __name__ == "__main__":
    main()
