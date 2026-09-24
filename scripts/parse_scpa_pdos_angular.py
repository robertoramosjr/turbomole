#!/usr/bin/env python3
"""
parse_scpa_pdos_angular.py

Converts a Multiwfn DOS_curve.txt export (SCPA-partitioned, angular-momentum
-resolved PDOS: s/p character per element) into a Veusz-importable dataset.

SCPA (Self-Consistent... method, see Multiwfn manual Sec. 2.4) was chosen
over Mulliken because it operates at the same basis-function granularity
(allowing angular-momentum-specific fragments, via the "l s"/"l p" commands
combined with atom-range "cond" restriction) while guaranteeing fragment
compositions never fall below 0% or exceed 100% -- unlike raw Mulliken,
which is the known source of the negative-PDOS artifact documented in
Sec. sec:dos of the Methods. Becke/Hirshfeld (used for the element-only
PDOS, Estagio 1b-alt2) cannot resolve angular momentum at all, since they
partition real space rather than basis functions (Multiwfn manual Sec. 2.4).

Usage:
    python ~/scripts/parse_scpa_pdos_angular.py --curve DOS_curve.txt \
        --output dos_scpa_angular_dataset.dat
"""

import argparse

# Column indices (0-based) in DOS_curve.txt, following the fragment
# definition order used for this molecule: 1=C-s, 2=C-p, 3=O-s, 4=O-p,
# 5=H-s, 6=H-p (fragments 7-10 undefined, always zero).
FRAGMENT_LABELS = ["pDOS_Cs", "pDOS_Cp", "pDOS_Os", "pDOS_Op", "pDOS_Hs", "pDOS_Hp"]


def parse_dos_curve(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 8:
                continue
            energy_ev = float(parts[0])
            tdos = float(parts[1])
            fragments = [float(parts[2 + i]) for i in range(6)]
            rows.append((energy_ev, tdos, fragments))
    if not rows:
        raise ValueError(f"No data rows parsed from {path}.")
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Convert a Multiwfn SCPA angular-momentum PDOS export "
                    "(6 fragments: C-s, C-p, O-s, O-p, H-s, H-p) to a Veusz dataset.")
    parser.add_argument("--curve", required=True, help="Path to DOS_curve.txt")
    parser.add_argument("--output", default="dos_scpa_angular_dataset.dat")
    args = parser.parse_args()

    rows = parse_dos_curve(args.curve)
    print(f"Parsed {len(rows)} points from {args.curve}")

    negative_counts = [0] * 6
    for _, _, fragments in rows:
        for i, value in enumerate(fragments):
            if value < 0:
                negative_counts[i] += 1

    for label, count in zip(FRAGMENT_LABELS, negative_counts):
        status = "OK, always non-negative" if count == 0 else f"WARNING: {count} negative points"
        print(f"  {label}: {status}")

    with open(args.output, "w") as f:
        header = "descriptor energy_eV,TDOS," + ",".join(FRAGMENT_LABELS)
        f.write(header + "\n")
        for energy_ev, tdos, fragments in rows:
            values = " ".join(f"{v:.6f}" for v in fragments)
            f.write(f"{energy_ev:.6f} {tdos:.6f} {values}\n")

    print(f"Veusz dataset saved -> {args.output}")


if __name__ == "__main__":
    main()