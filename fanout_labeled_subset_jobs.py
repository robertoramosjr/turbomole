#!/usr/bin/env python3
"""
fanout_labeled_subset_jobs.py

Given a production-level template folder (from
setup_hse06_tzvp_template.py) and the CSV written by
select_active_learning_subset.py, creates one job directory per
selected conformer: control/basis/auxbasis/mos copied from the
template, 'coord' copied from that conformer's r2SCAN-3c-optimized
geometry (Layer A). Same pattern as fanout_pbe0_jobs.py, but fans out
only the labeled subset instead of every Layer A survivor.

Usage:
    python ~/work_turbomole/scripts/fanout_labeled_subset_jobs.py \
        --template layer_b/hse06_template \
        --subset layer_b/initial_labeled_subset.csv \
        --r2scan3c-jobs layer_a/conformer_jobs \
        --output-dir layer_b/labeled_jobs
"""

import argparse
import csv
import os
import shutil


def main():
    parser = argparse.ArgumentParser(
        description="Fan out a production-level template into one job "
                    "directory per conformer in the active-learning "
                    "labeled subset.")
    parser.add_argument("--template", required=True)
    parser.add_argument("--subset", required=True,
                         help="select_active_learning_subset.py output CSV")
    parser.add_argument("--r2scan3c-jobs", required=True,
                         help="Parent folder of job_NNNN/coord dirs "
                              "(r2SCAN-3c-optimized geometries)")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    for fname in ("control", "basis", "mos"):
        if not os.path.isfile(os.path.join(args.template, fname)):
            raise SystemExit(f"Template '{args.template}' is missing '{fname}'.")
    # 'auxbasis' is optional here: Turbomole auto-generates it on the
    # first ridft/aoforce call when $rij is present but no matching
    # file exists yet (confirmed with the r2SCAN-3c template, which
    # never had one either) -- only RI-K templates (built via the
    # 'rijk'/'jkbas' define menu) pre-generate a real auxbasis file.
    optional_files = [f for f in ("auxbasis",)
                       if os.path.isfile(os.path.join(args.template, f))]

    with open(args.subset) as f:
        job_dirs = [row["job_dir"] for row in csv.DictReader(f)]
    if not job_dirs:
        raise SystemExit(f"No conformers listed in {args.subset}")

    if os.path.exists(args.output_dir):
        raise FileExistsError(f"{args.output_dir} already exists -- refusing to overwrite.")
    os.makedirs(args.output_dir)

    for jd in job_dirs:
        job_dir = os.path.join(args.output_dir, jd)
        os.makedirs(job_dir)
        for fname in ("control", "basis", "mos", *optional_files):
            shutil.copy(os.path.join(args.template, fname), os.path.join(job_dir, fname))
        shutil.copy(
            os.path.join(args.r2scan3c_jobs, jd, "coord"),
            os.path.join(job_dir, "coord"),
        )

    print(f"Created {len(job_dirs)} labeled-subset job directories under "
          f"{args.output_dir}/")


if __name__ == "__main__":
    main()
