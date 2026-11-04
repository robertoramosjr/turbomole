#!/usr/bin/env python3
"""Compara dois cálculos de frequência (ex.: aoforce × NumForce) pela Hessiana projetada.

Imprime frequências baixas, matriz de sobreposição |<a_i|b_j>| dos autovetores ponderados
por massa e cossenos principais dos subespaços dos k modos mais baixos.
Uso: compare_hessian_modes.py --a DIR_AOFORCE --b DIR_NUMFORCE [--nlow 8]
Cada DIR precisa de `coord` e `hessian`.
"""
import argparse
import os
import numpy as np
from tmhess import mass_weighted_modes

ap = argparse.ArgumentParser()
ap.add_argument('--a', required=True)
ap.add_argument('--b', required=True)
ap.add_argument('--nlow', type=int, default=8)
a = ap.parse_args()
fa, va = mass_weighted_modes(os.path.join(a.a, 'hessian'), os.path.join(a.a, 'coord'))
fb, vb = mass_weighted_modes(os.path.join(a.b, 'hessian'), os.path.join(a.b, 'coord'))
n = a.nlow
print('modos vibracionais:', len(fa), len(fb))
print('A:', np.round(fa[:n], 1))
print('B:', np.round(fb[:n], 1))
print('sobreposição |<B_i|A_j>| (linhas B, colunas A):')
print(np.round(np.abs(vb[:, :n].T @ va[:, :n]), 2))
for k in range(2, n):
    s = np.linalg.svd(vb[:, :k].T @ va[:, :k], compute_uv=False)
    print(f'subespaço dos {k} mais baixos, cossenos principais:', np.round(s, 2))
