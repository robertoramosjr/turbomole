#!/usr/bin/env python3
"""
extract_cosmo_energy.py

Extracts the COSMO(water) total energy (with outlying-charge
correction, the more accurate of the two values ridft prints) for a
set of conformers, from job_cosmo_ridft.sh's output.

Usage:
    python ~/work_turbomole/scripts/extract_cosmo_energy.py \
        --job-list job_dirs.txt \
        --cosmo-jobs layer_b/cosmo_jobs \
        --output layer_b/cosmo_energy.csv
"""

import argparse
import csv
import os
import re

OC_ENERGY_RE = re.compile(r"Total energy \+ OC corr\.\s*=\s*([-\d.]+)")


def main():
    parser = argparse.ArgumentParser(
        description="Extract COSMO total energy (with OC correction) from "
                    "ridft_cosmo.out files.")
    parser.add_argument("--job-list", required=True)
    parser.add_argument("--cosmo-jobs", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.job_list) as f:
        job_dirs = [line.strip() for line in f if line.strip()]

    rows = []
    for jd in job_dirs:
        path = os.path.join(args.cosmo_jobs, jd, "ridft_cosmo.out")
        if not os.path.isfile(path):
            print(f"WARNING: {path} not found, skipping")
            continue
        with open(path) as f:
            text = f.read()
        m = OC_ENERGY_RE.search(text)
        if not m:
            print(f"WARNING: no OC-corrected energy found in {path}")
            continue
        e = float(m.group(1))
        rows.append((jd, e))
        print(f"{jd}: E_cosmo(OC) = {e:.8f} Ha")

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["job_dir", "e_cosmo_oc_hartree"])
        for row in rows:
            writer.writerow([row[0], f"{row[1]:.8f}"])

    print(f"\nExtracted {len(rows)} / {len(job_dirs)} conformers.")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
