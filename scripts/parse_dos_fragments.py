#!/usr/bin/env python3
"""
parse_dos_fragments.py

Converts a Multiwfn DOS_curve.txt export (Becke/SCPA/Hirshfeld PDOS, any
number of fragments, energies in a.u.) into a Veusz-friendly dataset in
eV, with named fragment columns. Generalizes parse_becke_pdos.py (fixed
O/C/H) to arbitrary fragment definitions (e.g. ring_1, ring_2, azo_bridge
for azobenzene).

Multiwfn's DOS_curve.txt columns: energy | TDOS | frag.1 ... frag.10
(unused fragments are zero). The fragment columns do NOT sum to TDOS
(observed: sum/TDOS ~20.7 with ~+-2% variation across peaks, from the
Becke numerical integration), so by default each fragment is
renormalized: frag_k' = TDOS * frag_k / sum(frag). Use --no-normalize to
keep the raw columns.

Usage:
    python ~/scripts/parse_dos_fragments.py --curve DOS_curve.txt \
        --labels ring_1 ring_2 azo_bridge --output dos_becke_fragments.dat
"""

import argparse

import numpy as np

HARTREE_TO_EV = 27.211386245988


def main():
    parser = argparse.ArgumentParser(
        description="Convert a Multiwfn PDOS export (named fragments) to a Veusz dataset in eV.")
    parser.add_argument("--curve", required=True, help="Path to DOS_curve.txt")
    parser.add_argument("--labels", required=True, nargs="+",
                        help="Fragment names, in the order defined in Multiwfn (-1 menu)")
    parser.add_argument("--output", default="dos_becke_fragments.dat")
    parser.add_argument("--no-normalize", action="store_true",
                        help="Keep raw fragment columns instead of renormalizing to TDOS")
    args = parser.parse_args()

    data = np.loadtxt(args.curve)
    nfrag = len(args.labels)
    if data.shape[1] < 2 + nfrag:
        raise ValueError(f"{args.curve} has {data.shape[1]} columns, need >= {2 + nfrag}")

    energy_ev = data[:, 0] * HARTREE_TO_EV
    tdos = data[:, 1]
    frags = data[:, 2:2 + nfrag]

    if np.any(frags < 0) or np.any(tdos < 0):
        print("WARNING: negative values found -- check the Multiwfn export.")
    else:
        print("Confirmed: all TDOS and fragment values are non-negative.")

    total = frags.sum(axis=1)
    mask = tdos > 0.05 * tdos.max()
    ratio = total[mask] / tdos[mask]
    print(f"sum(frag)/TDOS where TDOS > 5% of max: min {ratio.min():.3f}  max {ratio.max():.3f} "
          f"(Multiwfn PDOS columns carry a constant ~20x scale vs TDOS plus ~+-2% Becke noise)")

    if not args.no_normalize:
        safe = np.where(total > 0, total, 1.0)
        frags = tdos[:, None] * frags / safe[:, None]
        print("Fragments renormalized so that sum(fragments) == TDOS at each energy.")

    with open(args.output, "w") as f:
        f.write("descriptor energy_eV,TDOS," + ",".join(f"pDOS_{l}" for l in args.labels) + "\n")
        for e, t, row in zip(energy_ev, tdos, frags):
            f.write(f"{e:.6f} {t:.6f} " + " ".join(f"{x:.6f}" for x in row) + "\n")
    print(f"Veusz dataset saved -> {args.output} ({len(energy_ev)} points)")


if __name__ == "__main__":
    main()
