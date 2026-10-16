#!/usr/bin/env python3
"""
batch_stda_s1.py

Runs sTDA (via std2, this project's compiled successor to the retired
grimme-lab/stda) on the r2SCAN-3c-converged wavefunction (Layer A,
already computed -- no new DFT) of every conformer in the Layer B
labeled subset, extracting the S1 excitation (energy, oscillator
strength) as B2's cheap proxy.

Pipeline per conformer: copy the r2SCAN-3c control/basis/coord/mos ->
`tm2molden` (answers its two interactive prompts with blank/default,
same as this project's earlier manual validation) -> `molden.input` ->
`std2 -f molden.input -ax 0.0 -e <threshold>` (ax=0.0 because r2SCAN is
a pure meta-GGA, 0% exact Hartree-Fock exchange -- NOT GFN2-xTB, which
was explicitly decided against as B2's proxy since it has no trustworthy
excited-state treatment; sTDA needs a real, if cheap, DFT reference,
which r2SCAN-3c already provides for free from Layer A) -> parse the
"state eV nm fL Rv(corr)" table for state 1 (S1, lowest energy).

Usage:
    python ~/work_turbomole/scripts/batch_stda_s1.py \
        --job-list job_dirs.txt \
        --r2scan3c-jobs layer_a/conformer_jobs \
        --energy-threshold 8.0 \
        --output layer_b/s1_stda.csv
"""

import argparse
import csv
import os
import re
import subprocess
import tempfile

TABLE_ROW_RE = re.compile(
    r"^\s*(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
)


def run_tm2molden(work_dir):
    result = subprocess.run(
        ["tm2molden"], cwd=work_dir, input="\n\n", capture_output=True,
        text=True, timeout=30,
    )
    return "ended normally" in result.stdout or "ended normally" in result.stderr


def run_std2(work_dir, energy_threshold, ax, std2_bin):
    result = subprocess.run(
        [std2_bin, "-f", "molden.input", "-ax", str(ax), "-e", str(energy_threshold)],
        cwd=work_dir, capture_output=True, text=True, timeout=120,
    )
    return result.stdout


def parse_s1(std2_stdout):
    lines = std2_stdout.splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("state") and "eV" in line and "fL" in line:
            header_idx = i
            break
    if header_idx is None:
        return None
    states = []
    for line in lines[header_idx + 1:]:
        m = TABLE_ROW_RE.match(line)
        if not m:
            break
        state, ev, nm, fl, rv = m.groups()
        states.append({"state": int(state), "energy_eV": float(ev), "osc": float(fl)})
    if not states:
        return None
    return min(states, key=lambda s: s["energy_eV"])


def main():
    parser = argparse.ArgumentParser(
        description="Batch sTDA S1 extraction on r2SCAN-3c wavefunctions.")
    parser.add_argument("--job-list", required=True)
    parser.add_argument("--r2scan3c-jobs", required=True,
                         help="Parent folder of job_NNNN/ with r2SCAN-3c "
                              "control+basis+coord+mos")
    parser.add_argument("--energy-threshold", type=float, default=8.0,
                         help="sTDA energy window in eV [default: 8.0]")
    parser.add_argument("--ax", type=float, default=0.0,
                         help="Fraction of exact Hartree-Fock exchange std2 "
                              "assumes for the reference wavefunction "
                              "[default: 0.0, appropriate for a pure "
                              "meta-GGA reference like r2SCAN-3c -- override "
                              "to test other values, e.g. 0.25 to match "
                              "HSE06's exact-exchange fraction]")
    parser.add_argument("--std2-bin", default=os.path.expanduser("~/local_bin/std2"))
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.job_list) as f:
        job_dirs = [line.strip() for line in f if line.strip()]

    results = []
    for jd in job_dirs:
        src = os.path.join(args.r2scan3c_jobs, jd)
        if not os.path.isdir(src):
            print(f"WARNING: {src} not found, skipping")
            continue

        with tempfile.TemporaryDirectory() as tmp:
            for fname in ("control", "basis", "coord", "mos"):
                src_file = os.path.join(src, fname)
                if os.path.isfile(src_file):
                    subprocess.run(["cp", src_file, tmp], check=True)

            if not run_tm2molden(tmp):
                print(f"WARNING: tm2molden failed for {jd}")
                continue
            if not os.path.isfile(os.path.join(tmp, "molden.input")):
                print(f"WARNING: no molden.input produced for {jd}")
                continue

            stdout = run_std2(tmp, args.energy_threshold, args.ax, args.std2_bin)
            s1 = parse_s1(stdout)
            if s1 is None:
                print(f"WARNING: could not parse sTDA output for {jd}")
                continue

            results.append((jd, s1["energy_eV"], s1["osc"]))
            print(f"{jd}: sTDA S1 = {s1['energy_eV']:.4f} eV, f = {s1['osc']:.5f}")

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["job_dir", "s1_energy_eV_stda", "s1_osc_stda"])
        for row in results:
            writer.writerow([row[0], f"{row[1]:.6f}", f"{row[2]:.6f}"])

    print(f"\nExtracted sTDA S1 for {len(results)} / {len(job_dirs)} conformers.")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
