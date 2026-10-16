#!/usr/bin/env python3
"""
select_b2_dft_confirmation_batch.py

Selects the batch of unconfirmed conformers to send to real
HSE06/TD-DFT (escf) for B2, for a protonation state where no per-
conformer delta-learning model survived calibration and every
non-ground-truth prediction is a flat offset applied to the sTDA
proxy (Sec. 9.20-9.21 of protonation_and_conformer_investigation.md
-- the neutral state, as of this project's B2 closure). Composition,
in priority order:

  1. --n-population-weighted conformers with the LARGEST population *
     f contribution (rank_b2_spectral_contribution.py's output) among
     those using the plain offset correction (method=surrogate_max_f,
     NOT uncertain_inflated) -- these are the conformers whose own
     confirmation individually moves the aggregate spectrum the most,
     and for which the offset is the ONLY information available (no
     per-conformer differentiation at all).
  2. --n-uncertain-inflated conformers flagged uncertain_inflated
     (intermediate sTDA selection-confidence zone, Sec. 9.19-9.20)
     with non-negligible population weight -- resolves state-selection
     ambiguity exactly where it matters most for the final spectrum.
  3. --n-diversity conformers chosen by farthest-point sampling
     (heavy-atom Kabsch RMSD, ignoring hydrogens) to maximize the
     minimum RMSD to everything already picked in steps 1-2 -- guards
     against the batch accidentally covering only one conformational
     family, regardless of population weight.

Always excludes: conformers marked 'excluded' (low sTDA selection
confidence, unresolved) and any job_dir in --already-resolved (manual
special-case overrides, e.g. job_0070/job_0012, which already have a
confirmed reassignment and do not need a new HSE06 calculation to
validate it).

Usage:
    python ~/work_turbomole/scripts/select_b2_dft_confirmation_batch.py \
        --contribution-rank layer_b/b2_spectral_contribution_rank.csv \
        --conformer-jobs layer_a/conformer_jobs \
        --already-resolved job_0070 job_0012 \
        --n-population-weighted 5 --n-uncertain-inflated 3 --n-diversity 2 \
        --output layer_b/b2_dft_confirmation_batch.csv
"""

import argparse
import csv
import os

import numpy as np


def read_coord_heavy(path):
    symbols, coords = [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 4 and parts[0].replace(".", "").replace("-", "").isdigit():
                x, y, z, sym = parts
                if sym.lower() != "h":
                    coords.append([float(x), float(y), float(z)])
                    symbols.append(sym.lower())
    return symbols, np.array(coords)


def kabsch_rmsd(p, q):
    p = p - p.mean(axis=0)
    q = q - q.mean(axis=0)
    h = p.T @ q
    u, s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    corr = np.diag([1, 1, d])
    r = vt.T @ corr @ u.T
    p_aligned = (r @ p.T).T
    return np.sqrt(np.mean(np.sum((p_aligned - q) ** 2, axis=1)))


def main():
    parser = argparse.ArgumentParser(
        description="Select a DFT-confirmation batch for a protonation "
                    "state whose B2 model is an offset only (no per-"
                    "conformer differentiation): top population-weighted "
                    "offset-only conformers + top uncertain-selection "
                    "conformers + forced geometric diversity.")
    parser.add_argument("--contribution-rank", required=True,
                         help="rank_b2_spectral_contribution.py output")
    parser.add_argument("--conformer-jobs", required=True,
                         help="Folder of job_NNNN/coord for RMSD diversity")
    parser.add_argument("--already-resolved", nargs="*", default=[],
                         help="job_dirs with a confirmed manual override -- "
                              "excluded, no new HSE06 needed")
    parser.add_argument("--n-population-weighted", type=int, default=5)
    parser.add_argument("--n-uncertain-inflated", type=int, default=3)
    parser.add_argument("--n-diversity", type=int, default=2)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    excluded_manual = set(args.already_resolved)

    with open(args.contribution_rank) as f:
        rows = list(csv.DictReader(f))
    pool = {r["job_dir"]: r for r in rows if r["job_dir"] not in excluded_manual}

    offset_only = sorted(
        (r for r in pool.values() if r["method"] == "surrogate_max_f"
         and r["uncertain_inflated"] != "True"),
        key=lambda r: float(r["contribution"]), reverse=True,
    )
    uncertain = sorted(
        (r for r in pool.values() if r["uncertain_inflated"] == "True"),
        key=lambda r: float(r["contribution"]), reverse=True,
    )

    selected = []
    reasons = {}
    for r in offset_only[:args.n_population_weighted]:
        selected.append(r["job_dir"])
        reasons[r["job_dir"]] = "top_population_weighted_offset_only"
    for r in uncertain[:args.n_uncertain_inflated]:
        if r["job_dir"] in selected:
            continue
        selected.append(r["job_dir"])
        reasons[r["job_dir"]] = "top_uncertain_inflated"

    print(f"Steps 1+2: {len(selected)} conformers selected by weight "
          f"({sum(1 for v in reasons.values() if v == 'top_population_weighted_offset_only')} "
          f"population-weighted offset-only + "
          f"{sum(1 for v in reasons.values() if v == 'top_uncertain_inflated')} "
          f"uncertain-inflated)")

    # Step 3: forced geometric diversity via FPS (heavy-atom Kabsch RMSD)
    coords = {}
    for jd in list(pool.keys()):
        path = os.path.join(args.conformer_jobs, jd, "coord")
        if os.path.isfile(path):
            _, xyz = read_coord_heavy(path)
            coords[jd] = xyz

    candidates = [jd for jd in pool if jd not in selected and jd in coords
                  and pool[jd]["method"] != "hse06_ground_truth"]
    for _ in range(args.n_diversity):
        if not candidates:
            break
        best_jd, best_min_rmsd = None, -1.0
        for jd in candidates:
            min_rmsd = min(kabsch_rmsd(coords[jd], coords[s]) for s in selected if s in coords)
            if min_rmsd > best_min_rmsd:
                best_min_rmsd = min_rmsd
                best_jd = jd
        selected.append(best_jd)
        reasons[best_jd] = f"forced_diversity_min_rmsd_{best_min_rmsd:.3f}A"
        candidates.remove(best_jd)
        print(f"Step 3: added {best_jd} (min heavy-atom RMSD to batch so far "
              f"= {best_min_rmsd:.3f} A)")

    print(f"\n=== Final batch (N={len(selected)}) ===")
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["job_dir", "population", "predicted_osc",
                           "contribution", "reason"])
        writer.writeheader()
        for jd in selected:
            r = pool[jd]
            writer.writerow({
                "job_dir": jd, "population": r["population"],
                "predicted_osc": r["predicted_osc"],
                "contribution": r["contribution"], "reason": reasons[jd],
            })
            print(f"  {jd}: pop={float(r['population']):.4f}, "
                  f"f={float(r['predicted_osc']):.4f}, "
                  f"contribution={float(r['contribution']):.5f} -- {reasons[jd]}")

    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
