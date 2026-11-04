"""Utilidades comuns: coord, $hessian, modos normais e massas (Turbomole 8.0.1).

Formato do `hessian`: linhas `<i> <bloco> v1..v5` (5 valores por linha, 3N*3N no total);
com NumForce o arquivo traz dois blocos (`$hessian` bruto e `$hessian (projected)`):
usa-se SEMPRE o projetado (último bloco).
"""
import numpy as np

BOHR_M = 0.529177210903e-10
AMU = 1.66053907e-27
EH = 4.3597447222e-18
C_CM = 2.99792458e10
MASS = {'h': 1.00794, 'c': 12.011, 'n': 14.0067, 'o': 15.9994, 's': 32.06}


def read_coord(path):
    rows = [l.split() for l in open(path).read().split('$coord')[1].split('$')[0].strip().splitlines()
            if len(l.split()) >= 4]
    return np.array([[float(x) for x in t[:3]] for t in rows]), [t[3].lower() for t in rows]


def masses3(symbols):
    return np.repeat([MASS[s] for s in symbols], 3)


def read_hessian(path, n3):
    vals = []
    for l in open(path):
        if l.startswith('$hessian (projected)'):
            vals = []
            continue
        if l.startswith('$hess'):
            continue
        if l.startswith('$end'):
            break
        vals.extend(float(x) for x in l.split()[2:])
    h = np.array(vals).reshape(n3, n3)
    return (h + h.T) / 2


def to_cm1(w):
    """autovalor de massa-ponderada (Eh/(amu bohr^2)) -> cm^-1 (negativo = imaginário)."""
    return np.sign(w) * np.sqrt(abs(w) * EH / (AMU * BOHR_M ** 2)) / (2 * np.pi * C_CM)


def cm1_to_k(f):
    """cm^-1 -> constante de força em Eh/(amu bohr^2)."""
    return np.sign(f) * (2 * np.pi * C_CM * f) ** 2 * AMU * BOHR_M ** 2 / EH


def mass_weighted_modes(hessian_path, coord_path, drop_below=1e-7):
    """(freq cm^-1, autovetores colunas) sem os 6 modos rígidos (|autovalor|<drop_below)."""
    X, sym = read_coord(coord_path)
    m = masses3(sym)
    h = read_hessian(hessian_path, 3 * len(sym))
    w, v = np.linalg.eigh(h / np.sqrt(np.outer(m, m)))
    keep = np.abs(w) > drop_below
    return to_cm1(w[keep]), v[:, keep]


def read_vibspectrum(path):
    import re
    out = []
    for l in open(path):
        if l[0] in '#$' or not l.strip():
            continue
        mt = re.search(r'(-?\d+\.\d+)\s+(\d+\.\d+)', l)
        if mt:
            out.append(float(mt.group(1)))
    return np.array(out)
