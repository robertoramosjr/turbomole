#!/usr/bin/env python3
"""
fanout_pbe0_jobs.py

Given a PBE0/def2-SVP RI-K + TD-DFT template folder (from
setup_pbe0_svp_template.py) and a parent folder of already-optimized
r2SCAN-3c conformer directories (layer_a/conformer_jobs/job_NNNN/, from
the earlier Layer A step), creates one new job directory per conformer:
control/basis/auxbasis/mos copied from the template, 'coord' copied
from that conformer's own final (jobex-optimized) geometry -- Turbomole
updates 'coord' in place during optimization, so job_NNNN/coord already
holds the r2SCAN-3c-optimized structure, no reconversion needed.

Usage:
    python ~/work_turbomole/scripts/fanout_pbe0_jobs.py \
        --template pbe0_template --r2scan3c-jobs layer_a/conformer_jobs \
        --output-dir layer_a/pbe0_jobs
"""

import argparse
import os
import shutil


def main():
    parser = argparse.ArgumentParser(
        description="Fan out a PBE0/def2-SVP template into one job "
                    "directory per r2SCAN-3c-optimized conformer.")
    parser.add_argument("--template", required=True,
                         help="Template dir with control/basis/auxbasis/mos "
                              "(from setup_pbe0_svp_template.py)")
    parser.add_argument("--r2scan3c-jobs", required=True,
                         help="Parent folder of job_NNNN/ dirs with "
                              "r2SCAN-3c-optimized 'coord' files")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    for fname in ("control", "basis", "auxbasis", "mos"):
        if not os.path.isfile(os.path.join(args.template, fname)):
            raise SystemExit(f"Template '{args.template}' is missing '{fname}'.")

    source_jobs = sorted(
        d for d in os.listdir(args.r2scan3c_jobs)
        if d.startswith("job_") and
        os.path.isfile(os.path.join(args.r2scan3c_jobs, d, "coord"))
    )
    if not source_jobs:
        raise SystemExit(f"No job_NNNN/coord found under {args.r2scan3c_jobs}")

    if os.path.exists(args.output_dir):
        raise FileExistsError(f"{args.output_dir} already exists -- refusing to overwrite.")
    os.makedirs(args.output_dir)

    for jd in source_jobs:
        job_dir = os.path.join(args.output_dir, jd)
        os.makedirs(job_dir)
        for fname in ("control", "basis", "auxbasis", "mos"):
            shutil.copy(os.path.join(args.template, fname), os.path.join(job_dir, fname))
        shutil.copy(
            os.path.join(args.r2scan3c_jobs, jd, "coord"),
            os.path.join(job_dir, "coord"),
        )

    print(f"Created {len(source_jobs)} PBE0 job directories under "
          f"{args.output_dir}/ (from r2SCAN-3c-optimized geometries in "
          f"{args.r2scan3c_jobs})")


if __name__ == "__main__":
    main()
