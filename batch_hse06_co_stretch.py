#!/usr/bin/env python3
"""
batch_hse06_co_stretch.py

Extracts the HSE06-D3(BJ)/def2-TZVP C=O stretch frequency (B1's
expensive target) for every conformer in the Layer B labeled subset,
via Turbomole's `aoforce2g98` (converts aoforce.out to the same
Gaussian-98-style format identify_co_stretch_mode.py already parses
for the GFN2-xTB proxy in batch_gfn2xtb_co_stretch.py) followed by the
same bond-vector-projection mode identification. Then joins the result
with that GFN2-xTB proxy table into one delta-learning-ready CSV
(target - proxy = the correction B1's GPR is trained to predict).

Usage:
    python ~/work_turbomole/scripts/batch_hse06_co_stretch.py \
        --subset layer_b/initial_labeled_subset.csv \
        --labeled-jobs layer_b/labeled_jobs \
        --r2scan3c-jobs layer_a/conformer_jobs \
        --gfn2xtb-table layer_b/gfn2xtb_co_stretch.csv \
        --output layer_b/b1_training_data.csv
"""

import argparse
import csv
import os
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
AOFORCE2G98 = os.path.expanduser(
    "~/software/TURBOMOLE/scripts/aoforce2g98")


def main():
    parser = argparse.ArgumentParser(
        description="Extract HSE06 C=O stretch frequencies for the labeled "
                    "subset and join with the GFN2-xTB proxy table.")
    parser.add_argument("--subset", required=True)
    parser.add_argument("--labeled-jobs", required=True,
                         help="Parent folder of job_NNNN/aoforce.out "
                              "(HSE06 production, from job_hse06_labeling*.sh)")
    parser.add_argument("--r2scan3c-jobs", required=True,
                         help="Parent folder of job_NNNN/structure.xyz "
                              "(r2SCAN-3c-optimized geometries)")
    parser.add_argument("--gfn2xtb-table", required=True,
                         help="batch_gfn2xtb_co_stretch.py output CSV")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.gfn2xtb_table) as f:
        gfn2xtb = {row["job_dir"]: float(row["nu_CO_gfn2xtb_cm1"])
                   for row in csv.DictReader(f)}

    with open(args.subset) as f:
        job_dirs = [row["job_dir"] for row in csv.DictReader(f)]

    results = []
    for jd in job_dirs:
        job_path = os.path.join(args.labeled_jobs, jd)
        aoforce_out = os.path.join(job_path, "aoforce.out")
        xyz_path = os.path.join(args.r2scan3c_jobs, jd, "structure.xyz")

        g98_result = subprocess.run(
            [AOFORCE2G98, aoforce_out], capture_output=True, text=True,
        )
        g98_path = os.path.join(job_path, "aoforce_g98.out")
        with open(g98_path, "w") as f:
            f.write(g98_result.stdout)

        identify_result = subprocess.run(
            ["python3", os.path.join(SCRIPT_DIR, "identify_co_stretch_mode.py"),
             "--g98", g98_path, "--xyz", xyz_path],
            capture_output=True, text=True,
        )
        freq = None
        for line in identify_result.stdout.splitlines():
            if line.startswith("C=O stretch mode:"):
                freq = float(line.split(",")[0].split(":")[1].strip().split()[0])

        if freq is None:
            print(f"WARNING: could not extract HSE06 C=O mode for {jd}")
            print(identify_result.stderr)
            continue
        if jd not in gfn2xtb:
            print(f"WARNING: {jd} has no GFN2-xTB proxy value, skipping")
            continue

        proxy = gfn2xtb[jd]
        results.append((jd, freq, proxy, freq - proxy))
        print(f"{jd}: HSE06 = {freq:.2f} cm-1, GFN2-xTB = {proxy:.2f} cm-1, "
              f"delta = {freq - proxy:+.2f} cm-1")

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["job_dir", "nu_CO_hse06_cm1", "nu_CO_gfn2xtb_cm1", "delta_cm1"])
        for row in results:
            writer.writerow(row)

    if results:
        deltas = [r[3] for r in results]
        mean_delta = sum(deltas) / len(deltas)
        spread = max(deltas) - min(deltas)
        print(f"\n{len(results)} / {len(job_dirs)} conformers matched.")
        print(f"Mean HSE06-GFN2xTB delta: {mean_delta:+.2f} cm-1, "
              f"spread across conformers: {spread:.2f} cm-1")
        print("(a small spread relative to the mean suggests the correction "
              "is close to a constant shift -- easy to delta-learn; a large "
              "spread means real conformer-to-conformer variation the GPR "
              "needs to capture, not just a systematic offset)")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
