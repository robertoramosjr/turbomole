#!/usr/bin/env python3
"""
batch_gfn2xtb_co_stretch.py

Runs GFN2-xTB Hessians (xtb --hess) on every conformer in the Layer B
labeled subset and extracts the C=O stretch frequency via
identify_co_stretch_mode.py's bond-projection method for each. This is
half of B1's delta-learning training data -- the cheap-proxy value,
matched against the HSE06/def2-TZVP value from the same conformer
(job_hse06_labeling.sh's aoforce output) once that finishes.

Usage:
    python ~/work_turbomole/scripts/batch_gfn2xtb_co_stretch.py \
        --subset layer_b/initial_labeled_subset.csv \
        --r2scan3c-jobs layer_a/conformer_jobs \
        --charge -1 \
        --output layer_b/gfn2xtb_co_stretch.csv
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    parser = argparse.ArgumentParser(
        description="Batch GFN2-xTB Hessian + C=O stretch identification "
                    "over a labeled-subset CSV.")
    parser.add_argument("--subset", required=True,
                         help="select_active_learning_subset.py output CSV")
    parser.add_argument("--r2scan3c-jobs", required=True,
                         help="Parent folder of job_NNNN/ dirs with "
                              "structure.xyz (r2SCAN-3c-optimized geometries)")
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--xtb-bin", default=os.path.expanduser(
        "~/.conda/envs/paper2_funnel/bin/xtb"))
    args = parser.parse_args()

    with open(args.subset) as f:
        job_dirs = [row["job_dir"] for row in csv.DictReader(f)]

    results = []
    for jd in job_dirs:
        xyz_path = os.path.join(args.r2scan3c_jobs, jd, "structure.xyz")
        if not os.path.isfile(xyz_path):
            print(f"WARNING: {xyz_path} not found, skipping {jd}", file=sys.stderr)
            continue

        with tempfile.TemporaryDirectory() as tmp:
            coord_path = os.path.join(tmp, "coord")
            subprocess.run(
                ["python3", os.path.join(SCRIPT_DIR, "xyz_to_coord.py"),
                 "--xyz", xyz_path, "--output", coord_path],
                check=True, capture_output=True,
            )
            hess_result = subprocess.run(
                [args.xtb_bin, "coord", "--gfn2", "--chrg", str(args.charge), "--hess"],
                cwd=tmp, capture_output=True, text=True,
            )
            g98_path = os.path.join(tmp, "g98.out")
            if hess_result.returncode != 0 or not os.path.isfile(g98_path):
                print(f"WARNING: xtb --hess failed for {jd}", file=sys.stderr)
                print(hess_result.stdout[-2000:], file=sys.stderr)
                continue

            identify_result = subprocess.run(
                ["python3", os.path.join(SCRIPT_DIR, "identify_co_stretch_mode.py"),
                 "--g98", g98_path, "--xyz", xyz_path],
                capture_output=True, text=True,
            )
            if identify_result.returncode != 0:
                print(f"WARNING: mode identification failed for {jd}", file=sys.stderr)
                print(identify_result.stderr, file=sys.stderr)
                continue

            freq = intensity = projection = None
            for line in identify_result.stdout.splitlines():
                if line.startswith("C=O stretch mode:"):
                    parts = line.split(",")
                    freq = float(parts[0].split(":")[1].strip().split()[0])
                    intensity = float(parts[1].split()[-1])
                    projection = float(parts[2].split()[-1])

            if freq is None:
                print(f"WARNING: could not parse mode for {jd}", file=sys.stderr)
                continue

            results.append((jd, freq, intensity, projection))
            print(f"{jd}: nu(C=O) = {freq:.2f} cm-1 (GFN2-xTB)")

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["job_dir", "nu_CO_gfn2xtb_cm1", "ir_intensity", "bond_projection"])
        for row in results:
            writer.writerow(row)

    print(f"\nProcessed {len(results)} / {len(job_dirs)} conformers -> {args.output}")


if __name__ == "__main__":
    main()
