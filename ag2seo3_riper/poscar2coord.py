#!/usr/bin/env python3
"""poscar2coord.py -- POSCAR -> TURBOMOLE periodic 'coord' (bohr) + cell parameters for riper.
Cell is rebuilt in the riper convention (a along x, b in the xy plane) from |a|,|b|,|c|,alpha,beta,gamma,
fractional coordinates are kept, so the structure is identical (rotated). Prints the $cell line in angstrom.
Usage: poscar2coord.py POSCAR OUTDIR"""
import sys
from pathlib import Path
import numpy as np
BOHR = 0.529177210903
L = Path(sys.argv[1]).read_text().splitlines()
s = float(L[1]); lat = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * s
el = L[5].split(); n = [int(x) for x in L[6].split()]
assert L[7].strip()[0] in "Dd", "expect Direct coordinates"
frac = np.array([[float(x) for x in l.split()[:3]] for l in L[8:8 + sum(n)]])
a, b, c = (np.linalg.norm(v) for v in lat)
al = np.degrees(np.arccos(lat[1] @ lat[2] / (b * c)))
be = np.degrees(np.arccos(lat[0] @ lat[2] / (a * c)))
ga = np.degrees(np.arccos(lat[0] @ lat[1] / (a * b)))
ca, cb, cg = (np.cos(np.radians(x)) for x in (al, be, ga)); sg = np.sin(np.radians(ga))
new = np.array([[a, 0, 0], [b * cg, b * sg, 0],
                [c * cb, c * (ca - cb * cg) / sg, c * np.sqrt(1 - cb**2 - ((ca - cb * cg) / sg)**2)]])
cart = frac @ new
out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
sym = [e.lower() for e, k in zip(el, n) for _ in range(k)]
with open(out / "coord", "w") as f:
    f.write("$coord\n")
    for (x, y, z), e in zip(cart / BOHR, sym): f.write(f"{x:22.14f}{y:22.14f}{z:22.14f}   {e}\n")
    f.write("$periodic 3\n$cell angs\n")
    f.write(f"  {a:.10f} {b:.10f} {c:.10f} {al:.8f} {be:.8f} {ga:.8f}\n$end\n")
print(f"a,b,c = {a:.5f} {b:.5f} {c:.5f}  alpha,beta,gamma = {al:.4f} {be:.4f} {ga:.4f}  V = {abs(np.linalg.det(new)):.3f} A^3  N = {len(sym)}")
