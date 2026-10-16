#!/usr/bin/env python3
"""
estimate_model_gap_impact.py

Quantifies, BY TYPE, what two known gaps of the artepilin C Layer B
cost in precision (protonation_and_conformer_investigation.md
Sec. 9.37):
  (1) monoanion B: the B1/B2 delta model is fitted on an O...H feature
      whose range across the 11 conformers is only ~0.0008 A;
  (2) monoanion A: no B2 (UV-Vis S1) model, real HSE06 for only the 20
      labeled conformers (26.6% of the population), B1 model Q2=0.28.

Types reported:
  T1  IR nu(C=O): spread of the real/predicted values, model vs. offset
      LOO error, and the effect on the population-weighted band.
  T2  UV-Vis: spread of S1 and of the FIRST BRIGHT transition across the
      real conformers, and the uncertainty of the population-weighted
      value from sampling only part of the ensemble.
  T3  pKa thermochemistry: sampling uncertainty of the per-state
      ensemble-average free energy (labeled subset vs. full ensemble),
      propagated to each deprotonation pair.
  T4  model validity: how much of the fitted "signal" a narrow-range
      feature can carry (target spread vs. feature range).

Sampling uncertainty for a state whose labeled subset covers only part
of the ensemble, with per-conformer values x_i (weighted sd sigma_w
over the labeled set) and full-ensemble weights w_i (sum = 1):
  random  : sigma_w * sqrt(sum_{unlabeled} w_i^2)   (unlabeled conformers
            scatter independently around the labeled mean)
  worst   : sigma_w * sum_{unlabeled} w_i           (ALL unlabeled shifted
            together by one sd -- fully correlated bound)
Both assume the labeled mean is unbiased for the ensemble; that is
tested separately with r2SCAN-3c energies, which exist for every
conformer (labeled-subset mean vs. full-ensemble mean).

Usage:
    python ~/work_turbomole/scripts/estimate_model_gap_impact.py \
        --project-root . --output model_gap_impact.dat
"""

import argparse
import csv
import math
import os
import re

import numpy as np

R_KCAL = 1.987204e-3
T_K = 298.15
HARTREE_TO_KCAL = 627.509474
KCAL_PER_PKA = R_KCAL * T_K * math.log(10)
ESCF_ROW_RE = re.compile(
    r"^\s*\|\s*(\d+)\s+a\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
    r"\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
)


def sdir(root, state):
    return os.path.join(root, f"artepillin_c_{state}_d3_hse06")


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def pop_full(root, state):
    rows = read_csv(os.path.join(sdir(root, state), "layer_a", "population_table_dft.csv"))
    p = {r["job_dir"]: float(r["r2scan3c_population"]) for r in rows}
    t = sum(p.values())
    return {k: v / t for k, v in p.items()}, rows


def wstats(x, w):
    w = np.asarray(w, float) / np.sum(w)
    x = np.asarray(x, float)
    m = float(np.sum(w * x))
    sd = float(math.sqrt(np.sum(w * (x - m) ** 2)))
    neff = float(1.0 / np.sum(w ** 2))
    return m, sd, neff


def sampling_unc(sd_labeled, w_full, labeled):
    unl = [w for j, w in w_full.items() if j not in labeled]
    return sd_labeled * math.sqrt(sum(w * w for w in unl)), sd_labeled * sum(unl)


def escf_states(path):
    st = []
    with open(path) as f:
        for line in f:
            m = ESCF_ROW_RE.match(line)
            if m:
                st.append((float(m.group(3)), float(m.group(7))))  # eV, f(length)
    return st


def first_bright(states, thr):
    for e, f in sorted(states):
        if f >= thr:
            return e, f
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Quantify by type the precision impact of the monoanion B "
                    "narrow-range model and the monoanion A missing B2 model.")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--bright-threshold", type=float, default=0.05)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    root = args.project_root
    out = []

    def say(s=""):
        print(s)
        out.append(s)

    # ------------------------------------------------------------ monoanion B
    say("################ MONOANION B (11/11 conformers = real HSE06) ################")
    wB, _ = pop_full(root, "monoanion_b")
    bdir = os.path.join(sdir(root, "monoanion_b"), "layer_b")
    oh = {r["job_dir"]: float(r["oh_distance_angstrom"]) for r in read_csv(os.path.join(bdir, "oh_distance_feature.csv"))}
    s1 = {r["job_dir"]: (float(r["s1_energy_eV"]), float(r["s1_osc_len"])) for r in read_csv(os.path.join(bdir, "s1_hse06.csv"))}
    b1 = {r["job_dir"]: float(r["nu_CO_hse06_cm1"]) for r in read_csv(os.path.join(bdir, "b1_training_data.csv"))}
    jds = sorted(oh)
    w = [wB[j] for j in jds]
    xs = np.array([oh[j] for j in jds])
    say(f"O...H feature range: {xs.min():.4f}-{xs.max():.4f} A (span {xs.max() - xs.min():.4f} A, N={len(xs)}); "
        f"distinct values: {sorted(set(np.round(xs, 4)))}")
    for name, vals, unit in (("E(S1)", [s1[j][0] for j in jds], "eV"),
                             ("f(S1)", [s1[j][1] for j in jds], ""),
                             ("nu(C=O)", [b1[j] for j in jds], "cm-1")):
        m, sd, neff = wstats(vals, w)
        say(f"  {name:8s}: weighted mean {m:.4f} {unit}, weighted sd {sd:.4f}, range {min(vals):.4f}..{max(vals):.4f} "
            f"(span {max(vals) - min(vals):.4f}), N_eff={neff:.1f}")
    Q2 = {"E": 0.5956, "logf": 0.4123}
    for key, name in (("E", "delta_energy_eV"), ("logf", "delta_log10_f")):
        tgt = [float(r[name]) for r in read_csv(os.path.join(bdir, "b2_osc_delta_target.csv"))]
        sd = float(np.std(tgt))
        loo_rmse = math.sqrt(1 - Q2[key]) * sd
        say(f"  model target {name}: sd {sd:.4f}; LOO RMSE of the O...H model ~ sqrt(1-Q2)*sd = {loo_rmse:.4f} "
            f"(Q2={Q2[key]}); slope*range would be extrapolating over {xs.max() - xs.min():.4f} A")
    # two-cluster check: is the feature just a family label?
    lo = [j for j in jds if oh[j] <= 2.2617]
    hi = [j for j in jds if oh[j] > 2.2617]
    fm = lambda js: float(np.mean([s1[j][1] for j in js]))
    say(f"  feature splits the set in two groups: {sorted(lo)} (O...H<=2.2617, mean f={fm(lo):.3f}) vs "
        f"{len(hi)} others (mean f={fm(hi):.3f}) -> the 'signal' is a group contrast, not a continuous relation")
    ratio = 10 ** (4.910 - 8.262)
    say(f"  population of B relative to A (tautomer ratio, calibrated pKa): [B]/[A] = {ratio:.2e} -> weight of B in the real "
        f"monoanion population ~{100 * ratio / (1 + ratio):.3f}%")
    e_sd = wstats([s1[j][0] for j in jds], w)[1]
    say(f"  => spectral error if B's conformer-to-conformer sd were wholly mis-modeled: {e_sd:.4f} eV in the B band; "
        f"diluted by the B weight in a monoanion mixture: {e_sd * ratio:.1e} eV. "
        f"NOTE: no model is applied to B (all 11 real), so the actual error from the model is 0.")

    # ------------------------------------------------------------ monoanion A
    say("\n################ MONOANION A (20/143 real HSE06 = 26.6% of population; B1 model on all 143) ################")
    wA, rowsA = pop_full(root, "monoanion_a")
    adir = os.path.join(sdir(root, "monoanion_a"), "layer_b")
    lab = sorted(os.listdir(os.path.join(adir, "labeled_jobs")))
    labset = set(lab)
    cov = sum(wA[j] for j in lab)
    say(f"labeled: {len(lab)}/{len(wA)}; population covered {cov:.3f}")
    # T1: B1
    tr = {r["job_dir"]: float(r["delta_cm1"]) for r in read_csv(os.path.join(adir, "b1_training_data.csv"))}
    nu = {r["job_dir"]: float(r["nu_CO_hse06_cm1"]) for r in read_csv(os.path.join(adir, "b1_training_data.csv"))}
    d = np.array([tr[j] for j in lab if j in tr])
    wl = [wA[j] for j in lab if j in tr]
    m, sd, neff = wstats(d, wl)
    say(f"T1 B1 delta(nu C=O) labeled: weighted mean {m:.2f} cm-1, weighted sd {sd:.2f}, range {d.min():.2f}..{d.max():.2f}")
    say(f"   LOO RMSE: offset-only ~ {np.std(d, ddof=1) * math.sqrt(len(d) / (len(d) - 1)):.2f} cm-1; Bader model (Q2=0.277) "
        f"~ {math.sqrt(1 - 0.277) * np.std(d):.2f} cm-1 -> the model buys {np.std(d, ddof=1) - math.sqrt(1 - 0.277) * np.std(d):.2f} cm-1 of accuracy")
    pred = {r["job_dir"]: float(r["predicted_delta_cm1"]) for r in read_csv(os.path.join(adir, "b1_predictions_ext.csv"))}
    allpred = dict(pred); allpred.update({j: tr[j] for j in tr})
    mf = sum(wA[j] * allpred[j] for j in wA if j in allpred) / sum(wA[j] for j in wA if j in allpred)
    say(f"   population-weighted delta over the FULL ensemble (real+model, {len([j for j in wA if j in allpred])} conformers): {mf:.2f} cm-1 "
        f"vs labeled-only {m:.2f} -> coverage extension changed the band position by {mf - m:+.2f} cm-1 (IR Lorentzian FWHM used: 10 cm-1)")
    # T2: UV-Vis from real escf
    thr = args.bright_threshold
    rowsUV = []
    for j in lab:
        st = escf_states(os.path.join(adir, "labeled_jobs", j, "escf.out"))
        fb = first_bright(st, thr)
        s1e, s1f = sorted(st)[0]
        rowsUV.append((j, s1e, s1f, fb))
    ok = [r for r in rowsUV if r[3] is not None]
    wS = [wA[r[0]] for r in rowsUV]
    m1, sd1, _ = wstats([r[1] for r in rowsUV], wS)
    say(f"T2 UV-Vis (real escf, N={len(rowsUV)}): S1 E weighted mean {m1:.3f} eV sd {sd1:.3f}; "
        f"f(S1) mean {np.mean([r[2] for r in rowsUV]):.5f} (dark: contributes ~nothing)")
    if ok:
        wB2 = [wA[r[0]] for r in ok]
        mb, sdb, nb = wstats([r[3][0] for r in ok], wB2)
        mf_, sdf, _ = wstats([r[3][1] for r in ok], wB2)
        say(f"   FIRST BRIGHT state (f>={thr}) found in {len(ok)}/{len(rowsUV)} conformers: E weighted mean {mb:.3f} eV, sd {sdb:.3f} "
            f"(range {min(r[3][0] for r in ok):.3f}..{max(r[3][0] for r in ok):.3f}); f mean {mf_:.3f}, sd {sdf:.3f}")
        rnd, wst = sampling_unc(sdb, wA, labset)
        say(f"   coverage uncertainty of that band position: random {rnd:.4f} eV, fully-correlated worst case {wst:.4f} eV")
        rndf, wstf = sampling_unc(sdf, wA, labset)
        say(f"   coverage uncertainty of its oscillator strength: random {rndf:.4f}, worst case {wstf:.4f}")
    # extra UV observables, less dependent on the brightness threshold
    strongest, centroid = [], []
    for j in lab:
        st = [x for x in escf_states(os.path.join(adir, "labeled_jobs", j, "escf.out")) if x[0] < 6.0]
        strongest.append(max(st, key=lambda x: x[1]))
        band = [x for x in st if 3.0 <= x[0] <= 5.5]
        centroid.append(sum(e * f for e, f in band) / sum(f for e, f in band))
    for name, vals in (("strongest transition <6 eV: E", [x[0] for x in strongest]),
                       ("strongest transition <6 eV: f", [x[1] for x in strongest]),
                       ("f-weighted centroid 3-5.5 eV", centroid)):
        mm, ss, _ = wstats(vals, [wA[j] for j in lab])
        rr, ww = sampling_unc(ss, wA, labset)
        say(f"   {name:32s}: weighted mean {mm:.3f}, sd {ss:.3f}, range {min(vals):.3f}..{max(vals):.3f}; "
            f"coverage unc random {rr:.4f} / worst {ww:.4f}")
    for thr2 in (0.02, 0.1):
        fb2 = [first_bright(escf_states(os.path.join(adir, "labeled_jobs", j, "escf.out")), thr2) for j in lab]
        v = [x[0] for x in fb2 if x]
        if len(v) == len(lab):
            mm, ss, _ = wstats(v, [wA[j] for j in lab])
            say(f"   first bright state, threshold f>={thr2}: E weighted mean {mm:.3f} eV, sd {ss:.3f}")
        else:
            say(f"   first bright state, threshold f>={thr2}: found in only {len(v)}/{len(lab)} conformers")
    # T3: thermochemistry for every state
    say("\n################ T3: pKa thermochemistry, sampling of the per-state ensemble average ################")
    unc = {}
    for state in ("neutral", "monoanion_a", "monoanion_b", "dianion"):
        sd_ = sdir(root, state)
        w_full, rows = pop_full(root, state)
        labs = set(os.listdir(os.path.join(sd_, "layer_b", "labeled_jobs")))
        freeh = {r["job_dir"]: float(r["g_thermal_correction_hartree"]) for r in read_csv(os.path.join(sd_, "layer_b", "freeh_thermal_correction.csv"))}
        cosmo = {r["job_dir"]: float(r["e_cosmo_oc_hartree"]) for r in read_csv(os.path.join(sd_, "layer_b", "cosmo_energy.csv"))}
        jl = sorted(set(freeh) & set(cosmo))
        g = np.array([(cosmo[j] + freeh[j]) * HARTREE_TO_KCAL for j in jl])
        wl = [w_full[j] for j in jl]
        mg, sdg, neffg = wstats(g, wl)
        rnd, wst = sampling_unc(sdg, w_full, set(jl))
        cov_ = sum(w_full[j] for j in jl)
        # representativeness with r2SCAN-3c energies (available for ALL conformers)
        e = {r["job_dir"]: float(r["r2scan3c_relative_energy_kcal"]) for r in rows}
        m_all = sum(w_full[j] * e[j] for j in w_full)
        m_lab = sum(w_full[j] * e[j] for j in jl) / sum(w_full[j] for j in jl)
        unc[state] = (rnd, wst, m_lab - m_all)
        say(f"{state:12s} labeled {len(jl):2d}/{len(w_full):3d} cov {cov_:.3f}; G_aq weighted sd {sdg:.3f} kcal/mol, N_eff {neffg:.1f}; "
            f"sampling unc random {rnd:.3f} / worst {wst:.3f} kcal/mol; r2SCAN-3c representativeness "
            f"(labeled mean - full mean) = {m_lab - m_all:+.3f} kcal/mol")
    say("\nPair propagation (pKa units; random = quadrature of the two states, worst = sum; "
        "representativeness = difference of the r2SCAN-3c subset-vs-full offsets):")
    for a, b, lab_ in (("neutral", "monoanion_a", "N->A"), ("neutral", "monoanion_b", "N->B"),
                       ("monoanion_a", "dianion", "A->D"), ("monoanion_b", "dianion", "B->D")):
        r = math.hypot(unc[a][0], unc[b][0]) / KCAL_PER_PKA
        wst = (unc[a][1] + unc[b][1]) / KCAL_PER_PKA
        rep = (unc[b][2] - unc[a][2]) / KCAL_PER_PKA
        say(f"  {lab_}: random {r:.3f}, worst-case {wst:.3f}, subset-representativeness {rep:+.3f} pKa")

    with open(args.output, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
