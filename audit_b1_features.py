#!/usr/bin/env python3
"""
audit_b1_features.py

Audits the B1 GPR's current feature vector before adding anything new
or attempting a pooled multi-state model: lists exactly which features
are in use (name, source parser, value range) and reports the simple
Pearson correlation of each feature IN ISOLATION against the delta
target, per protonation state -- this checks whether any individual
feature already carries linear signal that a GPR kernel (especially
one fit on N=11-20 points, as diagnose_b1_surrogate.py showed to be
unstable) might be diluting or, conversely, confirms none of them do.

Also explicitly checks whether bond-order (parse_bond_order.py, needs
pop.out from a Turbomole `dscf -proper`/`$pop` population run) or Bader
charge (parse_bader.py, needs ACF.dat from the external `bader` binary
on a `ridft -proper` density cube) data for the specific atom/bond of
the corrected C=O stretch mode is available for the labeled conformers
-- as opposed to only the two global geometric descriptors
(O...H distance, conjugation planarity) used so far. This project's
existing HSE06 labeling jobs (job_hse06_labeling*.sh) only ran
ridft+aoforce+escf, not a population-analysis step, so this is
expected to come back "not available yet" rather than "available but
unused" -- reported explicitly either way, not assumed.

Usage:
    python ~/work_turbomole/scripts/audit_b1_features.py \
        --training-data layer_b/b1_training_data.csv \
        --features layer_b/oh_distance_feature.csv layer_b/planarity_feature.csv \
        --labeled-jobs layer_b/labeled_jobs \
        --state-label neutral \
        --output layer_b/b1_feature_audit.dat
"""

import argparse
import csv
import os

import numpy as np
from scipy.stats import pearsonr


def load_data(training_path, feature_paths):
    with open(training_path) as f:
        training = {row["job_dir"]: float(row["delta_cm1"])
                    for row in csv.DictReader(f)}

    feature_values = {}
    feature_sources = {}
    for feat_path in feature_paths:
        with open(feat_path) as f:
            reader = csv.DictReader(f)
            cols = [c for c in reader.fieldnames if c != "job_dir"]
            for c in cols:
                feature_sources[c] = os.path.basename(feat_path)
            for row in reader:
                feature_values.setdefault(row["job_dir"], {})
                for c in cols:
                    feature_values[row["job_dir"]][c] = float(row[c])

    feature_names = list(feature_sources.keys())
    job_dirs = [jd for jd in training if jd in feature_values
                and all(c in feature_values[jd] for c in feature_names)]
    return job_dirs, training, feature_values, feature_names, feature_sources


def check_population_analysis_data(labeled_jobs_dir, job_dirs):
    has_pop = has_bader = 0
    for jd in job_dirs:
        job_path = os.path.join(labeled_jobs_dir, jd)
        if os.path.isfile(os.path.join(job_path, "pop.out")):
            has_pop += 1
        if os.path.isfile(os.path.join(job_path, "ACF.dat")):
            has_bader += 1
    return has_pop, has_bader, len(job_dirs)


def main():
    parser = argparse.ArgumentParser(
        description="Audit the B1 feature vector: sources, ranges, and "
                    "per-feature Pearson correlation with the delta target.")
    parser.add_argument("--training-data", required=True)
    parser.add_argument("--features", required=True, nargs="+")
    parser.add_argument("--labeled-jobs", required=True,
                         help="Parent folder of job_NNNN/ HSE06 production "
                              "dirs, to check for pop.out/ACF.dat availability")
    parser.add_argument("--state-label", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    job_dirs, training, feature_values, feature_names, feature_sources = \
        load_data(args.training_data, args.features)
    n = len(job_dirs)
    y = np.array([training[jd] for jd in job_dirs])

    print(f"=== {args.state_label} (N={n}) ===")
    print(f"Current feature vector ({len(feature_names)} features):")

    rows = []
    for name in feature_names:
        x = np.array([feature_values[jd][name] for jd in job_dirs])
        if n >= 3 and np.std(x) > 0:
            r, p = pearsonr(x, y)
        else:
            r, p = float("nan"), float("nan")
        print(f"  {name:28s} source={feature_sources[name]:28s} "
              f"range=[{x.min():.4f}, {x.max():.4f}]  "
              f"Pearson r={r:.4f} (p={p:.4f})")
        rows.append((name, feature_sources[name], x.min(), x.max(), r, p))

    has_pop, has_bader, total = check_population_analysis_data(args.labeled_jobs, job_dirs)
    print(f"\nBond-order data (parse_bond_order.py needs pop.out): "
          f"{has_pop} / {total} labeled conformers have it.")
    print(f"Bader charge data (parse_bader.py needs ACF.dat): "
          f"{has_bader} / {total} labeled conformers have it.")
    if has_pop == 0 and has_bader == 0:
        print("=> NEITHER is available for any labeled conformer yet -- "
              "the current feature vector is geometry-only (O...H distance, "
              "conjugation planarity), not because bond-order/Bader signal "
              "was tested and rejected, but because that data was never "
              "generated. job_hse06_labeling*.sh only ran ridft+aoforce+escf, "
              "no population-analysis step ($pop / dscf -proper / bader).")

    with open(args.output, "w") as f:
        f.write("descriptor state,feature_name,source,val_min,val_max,pearson_r,pearson_p\n")
        for name, src, vmin, vmax, r, p in rows:
            f.write(f"{args.state_label} {name} {src} {vmin:.4f} {vmax:.4f} "
                    f"{r:.4f} {p:.4f}\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
