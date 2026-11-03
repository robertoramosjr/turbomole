#!/usr/bin/env python3
"""Compare the CREST rotamer degeneracy g_i with the permutation images of equivalent atoms.

For an asymmetric conformer every element-preserving automorphism of the bond graph gives a distinct
labelled image, so the number of images is the SAME for all conformers.  If g_i were a symmetry/permutation
count it would be constant too.  This script (i) enumerates the automorphism group of the molecular graph
(bond graph from covalent radii, exhaustive backtracking), splitting it into heavy-atom images
(ring flip, gem-dimethyl swaps) and the full group (adds hydrogen permutations of the methyl groups), and
(ii) tabulates the distribution of g_i from water_ensemble/<state>/cre_members.

Usage:
  python degeneracy_vs_permutations.py --root ~/work_turbomole/artepillin_C_studies --out todo_v12/analysis/degeneracy_vs_permutations.csv
"""
import argparse
import collections
import csv
import os

import numpy as np

COV = {"h": 0.31, "c": 0.76, "o": 0.66, "n": 0.71}
STATES = {"N": "neutral", "A": "monoanion_a", "B": "monoanion_b", "D": "dianion"}


def read_coord(path):
    el, xyz, on = [], [], False
    for line in open(path):
        if line.startswith("$coord"):
            on = True
            continue
        if on and line.startswith("$"):
            break
        if on and line.split():
            p = line.split()
            xyz.append([float(v) * 0.529177 for v in p[:3]])
            el.append(p[3].lower())
    return el, np.array(xyz)


def graph(el, xyz):
    n = len(el)
    nb = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(xyz[i] - xyz[j]) < 1.25 * (COV[el[i]] + COV[el[j]]):
                nb[i].append(j)
                nb[j].append(i)
    return nb


def colours(el, nb):
    col = [hash(e) for e in el]
    for _ in range(len(el)):
        new = [hash((col[i], tuple(sorted(col[j] for j in nb[i])))) for i in range(len(el))]
        if len(set(new)) == len(set(col)):
            break
        col = new
    return col


def count_automorphisms(el, nb, only_heavy=False):
    """Number of graph automorphisms; if only_heavy, count distinct action on heavy atoms only."""
    n = len(el)
    col = colours(el, nb)
    order = []
    seen = set()
    start = max(range(n), key=lambda i: len(nb[i]))
    stack = [start]
    while stack:
        i = stack.pop(0)
        if i in seen:
            continue
        seen.add(i)
        order.append(i)
        stack.extend(j for j in nb[i] if j not in seen)
    images = set()
    cnt = [0]

    def rec(k, mp, used):
        if k == n:
            cnt[0] += 1
            if only_heavy:
                images.add(tuple(mp[i] for i in range(n) if el[i] != "h"))
            return
        i = order[k]
        for c in range(n):
            if c in used or col[c] != col[i]:
                continue
            ok = all((mp[j] in nb[c]) for j in nb[i] if j in mp) and \
                 all((j not in mp.values() or True) for j in nb[c])
            if ok and all(((c in nb[mp[j]]) == (j in nb[i])) for j in mp):
                mp[i] = c
                used.add(c)
                rec(k + 1, mp, used)
                del mp[i]
                used.discard(c)

    rec(0, {}, set())
    if only_heavy == "images":
        return sorted(images)
    return len(images) if only_heavy else cnt[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    root = os.path.expanduser(a.root)
    rows = []
    for st, d in STATES.items():
        el, xyz = read_coord(os.path.join(root, "water_ensemble", d, "coord"))
        nb = graph(el, xyz)
        n_full = count_automorphisms(el, nb)
        n_heavy = count_automorphisms(el, nb, only_heavy=True)
        g = [int(l.split()[0]) for l in list(open(os.path.join(root, "water_ensemble", d, "cre_members")))[1:] if l.strip()]
        c = collections.Counter(g)
        rows.append({"state": st, "n_struct": len(g), "aut_full": n_full, "images_heavy_only": n_heavy,
                     "g_min": min(g), "g_median": float(np.median(g)), "g_max": max(g),
                     "frac_g_eq_1": round(c[1] / len(g), 3), "n_distinct_g": len(c),
                     "sum_g": sum(g)})
        print(rows[-1])
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
