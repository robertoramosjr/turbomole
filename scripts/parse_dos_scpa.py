#!/usr/bin/env python3
"""
parse_scpa_pdos.py

Converts a Multiwfn DOS_curve.txt export (SCPA-partitioned PDOS, any
number of fragments) into a Veusz-importable dataset. Generalizes
parse_scpa_pdos_angular.py, which hardcoded exactly 6 fragments in a
fixed order (C-s, C-p, O-s, O-p, H-s, H-p) -- this version accepts
--labels for any fragment set, in the same order they were defined in
Multiwfn's "-1 Define fragments" menu.

SCPA was chosen over Mulliken because it operates at the same
basis-function granularity (allowing angular-momentum-specific
fragments, via the "l s"/"l p"/"l d" commands combined with atom-range
"cond" restriction, or a fragment covering the whole molecule) while
guaranteeing fragment compositions never fall below 0% or exceed
100% -- unlike raw Mulliken, which is the documented source of the
negative-PDOS artifact (Sec. sec:dos of the Methods), confirmed to
affect BOTH element-projected and whole-molecule angular-momentum
decompositions. Becke/Hirshfeld cannot resolve angular momentum at
all, since they partition real space rather than basis functions
(Multiwfn manual Sec. 2.4).

Usage:
    # element x angular-momentum PDOS (6 fragments: C-s, C-p, O-s, O-p, H-s, H-p)
    python ~/scripts/parse_scpa_pdos.py --curve DOS_curve.txt \
        --labels pDOS_Cs pDOS_Cp pDOS_Os pDOS_Op pDOS_Hs pDOS_Hp \
        --output dos_scpa_angular_dataset.dat

    # whole-molecule angular-momentum-only PDOS (3 fragments: s, p, d)
    python ~/scripts/parse_scpa_pdos.py --curve DOS_curve.txt \
        --labels DOS_s DOS_p DOS_d \
        --output dos_scpa_total_dataset.dat
"""

import argparse


def parse_dos_curve(path, n_fragments):
    rows = []
    with open(path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2 + n_fragments:
                continue
            energy_ev = float(parts[0])
            tdos = float(parts[1])
            fragments = [float(parts[2 + i]) for i in range(n_fragments)]
            rows.append((energy_ev, tdos, fragments))
    if not rows:
        raise ValueError(f"No data rows parsed from {path}.")
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Convert a Multiwfn SCPA PDOS export (any number of "
                    "fragments) to a Veusz dataset.")
    parser.add_argument("--curve", required=True, help="Path to DOS_curve.txt")
    parser.add_argument("--labels", required=True, nargs="+",
                         help="Descriptor column names, in the same order the "
                              "fragments were defined in Multiwfn (e.g. "
                              "pDOS_Cs pDOS_Cp pDOS_Os pDOS_Op pDOS_Hs pDOS_Hp, "
                              "or DOS_s DOS_p DOS_d)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    n_fragments = len(args.labels)
    rows = parse_dos_curve(args.curve, n_fragments)
    print(f"Parsed {len(rows)} points, {n_fragments} fragment(s) from {args.curve}")

    negative_counts = [0] * n_fragments
    for _, _, fragments in rows:
        for i, value in enumerate(fragments):
            if value < 0:
                negative_counts[i] += 1

    for label, count in zip(args.labels, negative_counts):
        status = "OK, always non-negative" if count == 0 else f"WARNING: {count} negative points"
        print(f"  {label}: {status}")

    with open(args.output, "w") as f:
        header = "descriptor energy_eV,TDOS," + ",".join(args.labels)
        f.write(header + "\n")
        for energy_ev, tdos, fragments in rows:
            values = " ".join(f"{v:.6f}" for v in fragments)
            f.write(f"{energy_ev:.6f} {tdos:.6f} {values}\n")

    print(f"Veusz dataset saved -> {args.output}")


if __name__ == "__main__":
    main()