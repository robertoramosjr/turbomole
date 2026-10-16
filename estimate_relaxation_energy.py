#!/usr/bin/env python3
"""
estimate_relaxation_energy.py

Newton-step estimate of the energy lost by NOT re-optimizing the
artepilin C conformers at HSE06 (protonation_and_conformer_investigation.md
Sec. 9.36). The conformers were optimized at r2SCAN-3c and only
single-point/Hessian-evaluated at HSE06-D3(BJ)/def2-TZVP, whereas the
pKa reference compounds were fully optimized at HSE06 -- so a
relaxation energy that differs between protonation states would bias
the deprotonation free energy in a way the reference offsets do NOT
absorb.

For each conformer: g = HSE06 Cartesian gradient at the frozen
geometry (rdgrad, job_rdgrad_relax_test.sh), H = the projected
Cartesian Hessian already computed there (aoforce, B1). Harmonic
relaxation energy over the vibrational (non-TR) subspace:

    dE = -1/2 * sum_k (v_k . g)^2 / lambda_k        (lambda_k > cutoff)

reported for a sweep of eigenvalue cutoffs (soft torsions/rotors have
small lambda and a huge Newton step -- the harmonic estimate is
unreliable there, so the sweep shows how much of dE comes from them),
plus the largest single-mode Newton step (bohr) as a validity check.

Per state, dE is population-averaged over the TESTED conformers
(weights = r2SCAN-3c Boltzmann populations, renormalized within the
tested set), and each deprotonation pair's shift is
    delta(deltaG) = <dE>_base - <dE>_acid ,
    delta(pKa)    = delta(deltaG) / (R T ln10).
Assumes the solvation energy is unchanged by the relaxation (gas-phase
gradient only) and that the tested top-population subset represents its
state.

Usage:
    python ~/work_turbomole/scripts/estimate_relaxation_energy.py \
        --joblist relax_test/joblist.txt --project-root . \
        --output relax_test/relaxation_energy.dat
"""

import argparse
import csv
import math
import os
import re

import numpy as np

HARTREE_TO_KCAL = 627.509474
R_KCAL = 1.987204e-3
T_K = 298.15
CUTOFFS = (1e-4, 5e-4, 1e-3, 5e-3)


def read_hessian(path, n3):
    h = np.zeros((n3, n3))
    with open(path) as f:
        for line in f:
            if line.startswith("$"):
                continue
            i = int(line[0:3]) - 1
            j = int(line[3:5]) - 1
            vals = [float(x) for x in line[5:].split()]
            h[i, j * 5: j * 5 + len(vals)] = vals
    return 0.5 * (h + h.T)


def read_last_gradient(path, natoms):
    with open(path) as f:
        lines = [l for l in f if l.strip()]
    starts = [i for i, l in enumerate(lines) if "cycle =" in l]
    s = starts[-1]
    grad_lines = lines[s + 1 + natoms: s + 1 + 2 * natoms]
    g = []
    for l in grad_lines:
        g.extend(float(x.replace("D", "E").replace("d", "e")) for x in l.split()[:3])
    return np.array(g)


def relaxation(h, g, cutoffs):
    w, v = np.linalg.eigh(h)
    proj = v.T @ g
    out = {}
    for c in cutoffs:
        mask = w > c
        out[c] = -0.5 * float(np.sum(proj[mask] ** 2 / w[mask])) * HARTREE_TO_KCAL
    mask = w > min(cutoffs)
    steps = np.abs(proj[mask] / w[mask])
    return out, float(steps.max()), int(np.sum(~mask)), w


def main():
    parser = argparse.ArgumentParser(
        description="Newton-step relaxation energy per conformer, its "
                    "population average per state, and the resulting pKa shifts.")
    parser.add_argument("--joblist", required=True,
                         help="relax_test joblist (<src>|<dst> per line, dst = "
                              "<state>_<job_dir>)")
    parser.add_argument("--project-root", required=True,
                         help="folder holding artepillin_c_<state>_d3_hse06/")
    parser.add_argument("--cutoff-for-average", type=float, default=1e-3,
                         help="Hessian-eigenvalue cutoff used for the per-state "
                              "average and pKa shift [default: 1e-3 Ha/bohr^2]")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = []
    with open(args.joblist) as f:
        for line in f:
            src, dst = line.strip().split("|")
            base = os.path.basename(dst)
            m = re.match(r"(.+)_(job_\d+)$", base)
            state, jd = m.group(1), m.group(2)
            gpath = os.path.join(dst, "gradient")
            if not os.path.isfile(gpath):
                print(f"WARNING: no gradient for {base}, skipping")
                continue
            with open(os.path.join(dst, "coord")) as cf:
                natoms = sum(1 for l in cf if len(l.split()) == 4 and l.split()[3].isalpha())
            g = read_last_gradient(gpath, natoms)
            h = read_hessian(os.path.join(dst, "hessian"), 3 * natoms)
            de, max_step, n_excl, w = relaxation(h, g, CUTOFFS)
            gmax = float(np.max(np.abs(g)))
            grms = float(np.sqrt(np.mean(g ** 2)))
            pop_path = os.path.join(args.project_root, f"artepillin_c_{state}_d3_hse06",
                                    "layer_a", "population_table_dft.csv")
            with open(pop_path) as pf:
                pop = {r["job_dir"]: float(r["r2scan3c_population"]) for r in csv.DictReader(pf)}[jd]
            rows.append({"state": state, "job_dir": jd, "pop": pop, "gmax": gmax,
                         "grms": grms, "max_step_bohr": max_step, "de": de,
                         "n_soft_excluded_at_min_cutoff": n_excl})

    print(f"{'state':12s} {'job':9s} {'gmax':>8s} {'grms':>8s} " +
          " ".join(f"dE(>{c:g})" .rjust(11) for c in CUTOFFS) + "  max_step(bohr)")
    for r in rows:
        print(f"{r['state']:12s} {r['job_dir']:9s} {r['gmax']:8.5f} {r['grms']:8.5f} " +
              " ".join(f"{r['de'][c]:11.3f}" for c in CUTOFFS) + f"  {r['max_step_bohr']:8.2f}")

    states = sorted({r["state"] for r in rows})
    avg = {}
    print(f"\n=== population-averaged relaxation energy (kcal/mol, cutoff {args.cutoff_for_average:g}) ===")
    for s in states:
        rs = [r for r in rows if r["state"] == s]
        wsum = sum(r["pop"] for r in rs)
        c = args.cutoff_for_average
        avg[s] = sum(r["pop"] * r["de"][c] for r in rs) / wsum
        spread = (min(r["de"][c] for r in rs), max(r["de"][c] for r in rs))
        print(f"{s:12s} N={len(rs)}  <dE> = {avg[s]:+.3f}   (per-conformer range {spread[0]:+.3f} .. {spread[1]:+.3f})")

    pairs = [("neutral", "monoanion_a", "N->A"), ("neutral", "monoanion_b", "N->B"),
             ("monoanion_a", "dianion", "A->D"), ("monoanion_b", "dianion", "B->D")]
    print("\n=== shift of the deprotonation free energy from NOT relaxing the artepilin geometries ===")
    print("(G_relaxed = G_frozen + dE with dE < 0, so pKa_relaxed = pKa_frozen + shift; "
          "add the shift to the RAW artepilin pKa, then apply the same reference offsets)")
    out_lines = []
    for a, b, label in pairs:
        if a in avg and b in avg:
            d = avg[b] - avg[a]
            dp = d / (R_KCAL * T_K * math.log(10))
            s = f"{label}: delta(deltaG) = {d:+.3f} kcal/mol = {dp:+.3f} pKa units"
            print(s); out_lines.append(s)

    with open(args.output, "w") as f:
        f.write("descriptor state,job_dir,pop,gmax,grms,max_newton_step_bohr," +
                ",".join(f"dE_kcal_cut{c:g}" for c in CUTOFFS) + "\n")
        for r in rows:
            f.write(f"{r['state']} {r['job_dir']} {r['pop']:.4f} {r['gmax']:.6f} {r['grms']:.6f} "
                    f"{r['max_step_bohr']:.3f} " + " ".join(f"{r['de'][c]:.4f}" for c in CUTOFFS) + "\n")
        f.write("\n# pair shifts (cutoff %g)\n" % args.cutoff_for_average)
        for s in out_lines:
            f.write("# " + s + "\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
