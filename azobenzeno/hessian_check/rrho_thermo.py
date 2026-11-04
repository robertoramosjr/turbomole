#!/usr/bin/env python3
"""Termodinâmica RRHO (gás ideal, 1 bar) a partir de vibspectrum + coord: ZPE, S, H, G.

Conferência do `freeh` (ZPE do Turbomole bate em 1e-7 Eh). |freq|<1 cm^-1 descartadas;
imaginárias entram como |freq| (avisa). sigma = nº de simetria rotacional (C1: 1; C2h: 2).
Uso: rrho_thermo.py --coord coord --vibspectrum vibspectrum --eel -572.67372535636 [--T 298.15] [--sigma 1]
"""
import argparse
import numpy as np
from tmhess import read_coord, read_vibspectrum, MASS, AMU, BOHR_M

kB, h, c, NA, Ha = 1.380649e-23, 6.62607015e-34, 2.99792458e10, 6.02214076e23, 4.3597447222e-18
KC = 627.5095
ap = argparse.ArgumentParser()
ap.add_argument('--coord', required=True)
ap.add_argument('--vibspectrum', required=True)
ap.add_argument('--eel', type=float, required=True)
ap.add_argument('--T', type=float, default=298.15)
ap.add_argument('--P', type=float, default=1e5, help='Pa')
ap.add_argument('--sigma', type=int, default=1)
a = ap.parse_args()
T = a.T
X, sym = read_coord(a.coord)
X = X * BOHR_M
m = np.array([MASS[s] for s in sym]) * AMU
Xc = X - (m[:, None] * X).sum(0) / m.sum()
I = sum(mi * (r @ r * np.eye(3) - np.outer(r, r)) for mi, r in zip(m, Xc))
Iev = np.linalg.eigvalsh(I)
f = read_vibspectrum(a.vibspectrum)
f = f[np.abs(f) > 1]
if (f < 0).any():
    print('AVISO: frequências imaginárias tratadas como |f|:', f[f < 0])
nu = np.abs(f) * c
x = h * nu / (kB * T)
zpe = 0.5 * h * nu.sum() / Ha
Evib = (h * nu * (0.5 + 1 / np.expm1(x))).sum() / Ha
Svib = (kB * (x / np.expm1(x) - np.log1p(-np.exp(-x)))).sum()
qt = (2 * np.pi * m.sum() * kB * T / h ** 2) ** 1.5 * kB * T / a.P
Str = kB * (np.log(qt) + 2.5)
qr = np.sqrt(np.pi) / a.sigma * np.prod(np.sqrt(8 * np.pi ** 2 * Iev * kB * T / h ** 2))
Srot = kB * (np.log(qr) + 1.5)
Eth = Evib + 3 * kB * T / Ha
Hh = a.eel + Eth + kB * T / Ha
S = Svib + Str + Srot
G = Hh - T * S / Ha
print(f'ZPE={zpe:.7f} Eh ({zpe*KC:.2f} kcal/mol)  S={S*NA/4.184:.2f} cal/mol/K  H-Eel={(Hh-a.eel)*KC:.2f}  '
      f'G-Eel={(G-a.eel)*KC:.2f} kcal/mol')
print(f'Eel={a.eel:.7f}  Eel+ZPE={a.eel+zpe:.7f}  H={Hh:.7f}  G={G:.7f} Eh')
