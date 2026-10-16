#!/usr/bin/env python3
"""
batch_extract_s1_hse06.py

Extracts the S1 (lowest-energy) TD-DFT excitation -- energy and
oscillator strength -- from each conformer's already-computed
HSE06/def2-TZVP escf.out (the SAME production labeling jobs used for
B1's Hessians, job_hse06_labeling*.sh already ran `escf` with 60
singlet states, matching Paper 1's window -- no new expensive DFT
needed for this half of B2's data). Reuses parse_excitations.py's
table-row regex for the "SUMMARY OF EXCITATION ENERGIES" block; S1 is
simply the lowest-energy row (state 1), since Turbomole lists TD-DFT
roots in ascending energy order.

Usage:
    python ~/work_turbomole/scripts/batch_extract_s1_hse06.py \
        --job-list job_dirs.txt \
        --labeled-jobs layer_b/labeled_jobs \
        --output layer_b/s1_hse06.csv
"""

import argparse
import csv
import os
import re

TABLE_ROW_RE = re.compile(
    r"^\s*\|\s*(\d+)\s+a\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
    r"\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
)


def extract_s1(escf_path):
    rows = []
    with open(escf_path) as f:
        for line in f:
            m = TABLE_ROW_RE.match(line)
            if m:
                state, e_eh, e_ev, e_cm1, e_nm, osc_vel, osc_len = m.groups()
                rows.append({
                    "state": int(state), "energy_eV": float(e_ev),
                    "osc_len": float(osc_len), "osc_vel": float(osc_vel),
                })
    if not rows:
        return None
    return min(rows, key=lambda r: r["energy_eV"])


def main():
    parser = argparse.ArgumentParser(
        description="Extract the S1 TD-DFT excitation (energy, oscillator "
                    "strength) from each conformer's HSE06 escf.out.")
    parser.add_argument("--job-list", required=True,
                         help="One job_dir name per line (matching --labeled-jobs "
                              "subdirectory names)")
    parser.add_argument("--labeled-jobs", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.job_list) as f:
        job_dirs = [line.strip() for line in f if line.strip()]

    results = []
    for jd in job_dirs:
        escf_path = os.path.join(args.labeled_jobs, jd, "escf.out")
        if not os.path.isfile(escf_path):
            print(f"WARNING: {escf_path} not found, skipping")
            continue
        s1 = extract_s1(escf_path)
        if s1 is None:
            print(f"WARNING: no excitation table found in {escf_path}, skipping")
            continue
        results.append((jd, s1["energy_eV"], s1["osc_len"], s1["osc_vel"]))
        print(f"{jd}: S1 = {s1['energy_eV']:.4f} eV, f(len) = {s1['osc_len']:.5f}")

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["job_dir", "s1_energy_eV", "s1_osc_len", "s1_osc_vel"])
        for row in results:
            writer.writerow([row[0], f"{row[1]:.6f}", f"{row[2]:.6f}", f"{row[3]:.6f}"])

    print(f"\nExtracted S1 data for {len(results)} / {len(job_dirs)} conformers.")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
