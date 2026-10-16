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


def parse_dos_curve(path, n_fragments, opdos=False, energy_au=False, scale=1.0):
    """Fixed-width parse (Multiwfn glues columns when a value is negative,
    e.g. ``29.468172-1019.275168``, so split() is unsafe). Layout: energy
    (10 chars) + 12-char fields. With >=2 fragments defined, column 3 is the
    OPDOS between fragments 1 and 2 (negative by nature, NOT a fragment):
    pass opdos=True to skip it."""
    HARTREE_TO_EV = 27.211386245988
    rows = []
    skip = 1 if opdos else 0
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if len(line) < 22:
                continue
            vals = [float(line[10 + i:22 + i]) for i in range(0, len(line) - 10, 12)
                    if line[10 + i:22 + i].strip()]
            if len(vals) < 1 + skip + n_fragments:
                continue
            energy = float(line[:10]) * (HARTREE_TO_EV if energy_au else 1.0)
            tdos = vals[0] * scale
            fragments = [v * scale for v in vals[1 + skip:1 + skip + n_fragments]]
            rows.append((energy, tdos, fragments))
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
    parser.add_argument("--opdos", action="store_true",
                        help="DOS_curve.txt has an OPDOS column after TDOS "
                             "(any run with >=2 fragments defined): skip it.")
    parser.add_argument("--energy-au", action="store_true",
                        help="Energy column is in a.u. (convert to eV).")
    parser.add_argument("--scale", type=float, default=1.0,
                        help="Multiply DOS columns (10 matches the Becke "
                             "parser convention, scale ratio 0.1).")
    args = parser.parse_args()

    n_fragments = len(args.labels)
    rows = parse_dos_curve(args.curve, n_fragments, args.opdos,
                           args.energy_au, args.scale)
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