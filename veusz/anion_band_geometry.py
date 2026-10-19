"""H3: conjugation of the carboxyl(ate) with the vinyl/ring in vacuum vs water-reoptimised geometries.

Atoms are identified by graph (covalent-radius bonds), not by index: carboxyl C = C with 2 O and 1 C
neighbour; C_alpha = that C neighbour; C_beta = other C neighbour of C_alpha; C_ipso = other C neighbour
of C_beta; C_ortho = aromatic C neighbour of C_ipso.  Vacuum: HSE06-labelled conformers (weights =
r2scan3c_population); water: every reopt dir with GEO_OPT_CONVERGED (unweighted; cNNNN numbering is
not the vacuum job numbering, so only distributions are compared).
Run with: LC_ALL=C ~/.conda/envs/veusz_env/bin/python anion_band_geometry.py
"""
import argparse, csv, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from anion_band_coverage import SDIR, hse_jobs, read_pop  # noqa: E402

BOHR = 0.529177
RAD = {"h": 0.31, "c": 0.76, "o": 0.66}


def read_coord(path):
    el, xyz = [], []
    for ln in open(path):
        t = ln.split()
        if len(t) == 4 and t[3].isalpha():
            xyz.append([float(x) * BOHR for x in t[:3]]); el.append(t[3].lower())
        elif ln.startswith("$") and xyz:
            break
    return el, np.array(xyz)


def neighbours(el, xyz):
    n = len(el)
    nb = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(xyz[i] - xyz[j]) < 1.25 * (RAD[el[i]] + RAD[el[j]]):
                nb[i].append(j); nb[j].append(i)
    return nb


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def twist(a):  # deviation from planarity, 0..90 deg (cis/trans equivalent)
    a = abs(a) % 180
    return min(a, 180 - a)


def measure(path):
    el, x = read_coord(path)
    nb = neighbours(el, x)
    C = lambda i: [j for j in nb[i] if el[j] == "c"]
    O = lambda i: [j for j in nb[i] if el[j] == "o"]
    ccarb = [i for i in range(len(el)) if el[i] == "c" and len(O(i)) == 2 and len(C(i)) == 1]
    if len(ccarb) != 1:
        raise ValueError(f"{path}: {len(ccarb)} carboxyl carbons")
    c0 = ccarb[0]; ca = C(c0)[0]
    cb = [j for j in C(ca) if j != c0]
    if len(cb) != 1:
        raise ValueError(f"{path}: C_alpha has {len(cb)} C neighbours")
    cb = cb[0]; ci = [j for j in C(cb) if j != ca][0]
    co = [j for j in C(ci) if len([k for k in nb[j] if el[k] == "h"]) <= 1]
    d = np.linalg.norm
    o1, o2 = O(c0)
    return dict(tw_OCCa_Cb=min(twist(dihedral(x[o], x[c0], x[ca], x[cb])) for o in (o1, o2)),
                tw_Ca_Cb_Ci_Co=min(twist(dihedral(x[ca], x[cb], x[ci], x[j])) for j in co),
                tw_Cc_Ca_Cb_Ci=twist(dihedral(x[c0], x[ca], x[cb], x[ci])),
                r_CaCb=d(x[ca] - x[cb]), r_CcCa=d(x[c0] - x[ca]), r_CbCi=d(x[cb] - x[ci]),
                r_CO_diff=abs(d(x[c0] - x[o1]) - d(x[c0] - x[o2])))


def vacuum_paths(d):
    pop = read_pop(d)
    for j in hse_jobs(d) & set(pop):
        for sb in ("labeled_jobs", "confirmation_batch_jobs", "random_control_batch_jobs"):
            p = os.path.join(d, "layer_b", sb, j, "coord")
            if os.path.isfile(p):
                yield j, pop[j], p; break


def water_paths(root, state):
    base = os.path.join(root, "reopt", state)
    for c in sorted(os.listdir(base)):
        if os.path.isfile(os.path.join(base, c, "GEO_OPT_CONVERGED")):
            yield c, 1.0, os.path.join(base, c, "coord")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--root", default=os.path.abspath(os.path.join(here, "..", "..")))
    ap.add_argument("--out", default=os.path.join(here, "..", "Figures", "data", "anion_band_geometry.csv"))
    a = ap.parse_args()
    rows = []
    for s, name in SDIR.items():
        d = os.path.join(a.root, f"artepillin_c_{name}_d3_hse06")
        for med, it in (("vacuum", vacuum_paths(d)), ("water_reopt", water_paths(a.root, name))):
            for cid, w, p in it:
                rows.append(dict(state=s, medium=med, conf=cid, weight=w, **measure(p)))
    keys = list(rows[0])
    with open(a.out, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys); wr.writeheader()
        wr.writerows({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()} for r in rows)
    for s in SDIR:
        for med in ("vacuum", "water_reopt"):
            sel = [r for r in rows if r["state"] == s and r["medium"] == med]
            w = np.array([r["weight"] for r in sel]); w = w / w.sum()
            print(f"{s} {med:11s} n={len(sel):3d} " + " ".join(
                f"{k}={np.sum(w * [r[k] for r in sel]):.2f}" for k in keys[4:]))


if __name__ == "__main__":
    main()
