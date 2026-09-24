#!/usr/bin/env python3
"""
compare_functionals.py

Extracts the dominant TD-DFT excited state (highest oscillator
strength) from multiple escf.out files (one per functional) and
builds a comparison table -- the "compare HSE06 with a global hybrid
and a range-separated hybrid for the key excitation energies" part of
the robustness requirement, with a fixed conformer/geometry.

Usage:
    python ~/scripts/compare_functionals.py \
        --run hse06-d3bj:artepillin_C_d3/escf_d3.out \
        --run b3lyp:artepillin_C_d3_b3lyp/escf_b3lyp.out \
        --run cam-b3lyp:artepillin_C_d3_camb3lyp/escf_camb3lyp.out \
        --output functional_comparison.csv
"""

import argparse
import re

TABLE_ROW_RE = re.compile(
    r"^\s*\|\s*(\d+)\s+a\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
)


def parse_escf_table(escf_path):
    rows = []
    with open(escf_path, "r") as f:
        for line in f:
            match = TABLE_ROW_RE.match(line)
            if match:
                state, e_eh, e_ev, e_cm1, e_nm, osc_vel, osc_len = match.groups()
                rows.append({
                    "state": int(state), "energy_eV": float(e_ev),
                    "energy_nm": float(e_nm), "osc_len": float(osc_len),
                })
    if not rows:
        raise ValueError(f"No excitation rows found in {escf_path}")
    return rows


def main():
    parser = argparse.ArgumentParser(description="Compare dominant TD-DFT excitation across functionals.")
    parser.add_argument("--run", action="append", required=True,
                         metavar="LABEL:PATH",
                         help="functional label and path to its escf.out, e.g. b3lyp:path/escf.out. Repeat.")
    parser.add_argument("--output", default="functional_comparison.csv")
    args = parser.parse_args()

    results = []
    for entry in args.run:
        label, path = entry.split(":", 1)
        rows = parse_escf_table(path)
        dominant = max(rows, key=lambda r: r["osc_len"])
        results.append({"functional": label, **dominant})
        print(f"[{label}] dominant state {dominant['state']}: "
              f"{dominant['energy_eV']:.4f} eV ({dominant['energy_nm']:.1f} nm), "
              f"f = {dominant['osc_len']:.5f}")

    with open(args.output, "w") as f:
        f.write("functional,state,energy_eV,energy_nm,osc_len\n")
        for r in results:
            f.write(f"{r['functional']},{r['state']},{r['energy_eV']:.4f},"
                     f"{r['energy_nm']:.2f},{r['osc_len']:.5f}\n")
    print(f"\nComparison table saved -> {args.output}")


if __name__ == "__main__":
    main()