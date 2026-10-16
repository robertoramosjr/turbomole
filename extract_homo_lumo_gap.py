#!/usr/bin/env python3
"""
extract_homo_lumo_gap.py

Extracts the HSE06/def2-TZVP HOMO-LUMO gap (already printed by
`ridft` in every B1/B2 production labeling job -- no new calculation
needed) for a list of conformers, as one candidate descriptor for pKa
calibration (protonation_and_conformer_investigation.md Sec. 9.30).

Usage:
    python ~/work_turbomole/scripts/extract_homo_lumo_gap.py \
        --job-list job_dirs.txt \
        --labeled-jobs layer_b/labeled_jobs \
        --output layer_b/homo_lumo_gap.csv
"""

import argparse
import csv
import os
import re

GAP_RE = re.compile(r"HOMO-LUMO gap:\s*[-+\d.]+\s*H\s*=\s*([-+\d.]+)\s*eV")
HOMO_RE = re.compile(r"HOMO\s*:\s*[-+\d.]+\s*H\s*=\s*([-+\d.]+)\s*eV")
LUMO_RE = re.compile(r"LUMO\s*:\s*[-+\d.]+\s*H\s*=\s*([-+\d.]+)\s*eV")


def extract(ridft_path):
    with open(ridft_path) as f:
        text = f.read()
    gap = GAP_RE.search(text)
    homo = HOMO_RE.search(text)
    lumo = LUMO_RE.search(text)
    if not (gap and homo and lumo):
        return None
    return float(homo.group(1)), float(lumo.group(1)), float(gap.group(1))


def main():
    parser = argparse.ArgumentParser(
        description="Extract HOMO-LUMO gap (eV) from existing ridft.out "
                    "files -- no new DFT needed.")
    parser.add_argument("--job-list", required=True)
    parser.add_argument("--labeled-jobs", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.job_list) as f:
        job_dirs = [line.strip() for line in f if line.strip()]

    rows = []
    for jd in job_dirs:
        ridft_path = os.path.join(args.labeled_jobs, jd, "ridft.out")
        if not os.path.isfile(ridft_path):
            print(f"WARNING: {ridft_path} not found, skipping")
            continue
        result = extract(ridft_path)
        if result is None:
            print(f"WARNING: could not find HOMO-LUMO block in {ridft_path}")
            continue
        homo, lumo, gap = result
        rows.append((jd, homo, lumo, gap))
        print(f"{jd}: HOMO={homo:.4f} eV, LUMO={lumo:.4f} eV, gap={gap:.4f} eV")

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["job_dir", "homo_eV", "lumo_eV", "homo_lumo_gap_eV"])
        for row in rows:
            writer.writerow([row[0], f"{row[1]:.6f}", f"{row[2]:.6f}", f"{row[3]:.6f}"])

    print(f"\nExtracted {len(rows)} / {len(job_dirs)} conformers.")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
