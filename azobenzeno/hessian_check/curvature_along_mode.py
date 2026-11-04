#!/usr/bin/env python3
"""Teste de curvatura ao longo de UM modo normal: decide se uma imaginária do aoforce é real.

  gen:     cria DIR/<tag>_{p,m}<q>/ com coord deslocada só ao longo do modo
           (dx = q*v/sqrt(m), q em amu^1/2*bohr), copiando control/basis/auxbasis/mos de --base
           (--base deve ser o ponto estacionário já convergido, com a MESMA grade do teste) e
           escrevendo DIR/cases.txt (um caso por linha, p/ job array). Rodar em cada caso:
           actual -r; ridft; rdgrad   (rdgrad grava o gradiente em $grad do `control`).
  analyze: lê ridft.out e $grad e imprime ΔE(±q), dE/dq e ω = sqrt(k) com k=(g+ - g-)/2q,
           mais a extrapolação Richardson q->0 (k(q)=k0+a q^2) com os dois menores q.
Modo escolhido: --mode-index i (0 = mais baixo, contando só os modos não rígidos) do --hessian.
Uso:
  curvature_along_mode.py gen --base DIR0 --hessian H --coord C --mode-index 0 --tag g7_ao \
        --q 1 2 3 --outdir curvA [--gridsize 7]
  curvature_along_mode.py analyze --outdir curvA --tag g7_ao --q 1 2 3 --hessian H --coord C --mode-index 0
"""
import argparse
import os
import re
import shutil
import numpy as np
from tmhess import mass_weighted_modes, read_coord, masses3, cm1_to_k, to_cm1


def mode_vector(a):
    _, v = mass_weighted_modes(a.hessian, a.coord)
    return v[:, a.mode_index]


def gen(a):
    X0, sym = read_coord(a.coord)
    m = masses3(sym)
    v = mode_vector(a)
    txt = open(os.path.join(a.base, 'coord')).read()
    head, rest = txt.split('$coord')[0], txt.split('$coord')[1]
    tail = '$' + '$'.join(rest.split('$')[1:])
    cases = []
    for q in a.q:
        for s, sn in ((1, 'p'), (-1, 'm')):
            name = f'{a.tag}_{sn}{q:g}'
            d = os.path.join(a.outdir, name)
            os.makedirs(d, exist_ok=True)
            for f in ('auxbasis', 'basis', 'control', 'mos'):
                shutil.copy(os.path.join(a.base, f), d)
            if a.gridsize:
                c = open(os.path.join(d, 'control')).read()
                c = re.sub(r'(gridsize\s+)\S+', rf'\g<1>{a.gridsize}', c)
                open(os.path.join(d, 'control'), 'w').write(c)
            X = X0 + (s * q * v / np.sqrt(m)).reshape(-1, 3)
            lines = ['%20.14f%20.14f%20.14f      %s' % (*X[i], sym[i]) for i in range(len(sym))]
            open(os.path.join(d, 'coord'), 'w').write(head + '$coord\n' + '\n'.join(lines) + '\n' + tail)
            cases.append(name)
    with open(os.path.join(a.outdir, 'cases.txt'), 'a') as f:
        f.write('\n'.join(cases) + '\n')
    print(len(cases), 'casos em', a.outdir)


def energy(d):
    return float(re.findall(r'total energy\s+=\s+(-?\d+\.\d+)', open(os.path.join(d, 'ridft.out')).read())[-1])


def grad(d, nat):
    t = open(os.path.join(d, 'control')).read().split('$grad')[1].split('$')[0].splitlines()
    cyc = [i for i, l in enumerate(t) if 'cycle' in l][-1]
    rows = [l for l in t[cyc + 1:] if l.strip()]
    return np.array([[float(x.replace('D', 'E')) for x in l.split()[:3]] for l in rows[nat:2 * nat]]).ravel()


def analyze(a):
    X0, sym = read_coord(a.coord)
    m = masses3(sym)
    u = mode_vector(a) / np.sqrt(m)
    ks = {}
    for q in a.q:
        dp = os.path.join(a.outdir, f'{a.tag}_p{q:g}')
        dm = os.path.join(a.outdir, f'{a.tag}_m{q:g}')
        gp, gm = grad(dp, len(sym)) @ u, grad(dm, len(sym)) @ u
        ks[q] = (gp - gm) / (2 * q)
        print(f'q=±{q:g}: E+={energy(dp):.8f} E-={energy(dm):.8f}  dE/dq(+)={gp:+.3e} dE/dq(-)={gm:+.3e}'
              f'  k={ks[q]:+.3e} Eh/(amu bohr2)  w={to_cm1(ks[q]):+.1f} cm-1')
    q1, q2 = sorted(a.q)[:2]
    aq = (ks[q2] - ks[q1]) / (q2 ** 2 - q1 ** 2)
    k0 = ks[q1] - aq * q1 ** 2
    print(f'extrapolação q->0 (k=k0+a q^2): k0={k0:+.3e}  w0={to_cm1(k0):+.1f} cm-1')


ap = argparse.ArgumentParser()
sub = ap.add_subparsers(dest='cmd', required=True)
for name in ('gen', 'analyze'):
    p = sub.add_parser(name)
    p.add_argument('--hessian', required=True)
    p.add_argument('--coord', required=True)
    p.add_argument('--mode-index', type=int, default=0, dest='mode_index')
    p.add_argument('--tag', required=True)
    p.add_argument('--q', type=float, nargs='+', default=[1.0, 2.0, 3.0])
    p.add_argument('--outdir', required=True)
    if name == 'gen':
        p.add_argument('--base', required=True)
        p.add_argument('--gridsize', default=None)
a = ap.parse_args()
{'gen': gen, 'analyze': analyze}[a.cmd](a)
