#!/usr/bin/env python3
"""
parse_homo_lumo.py

Parses the "HOMO-LUMO Separation" block printed by Turbomole's -proper
runs (ridft -proper / dscf -proper) into a single-row Veusz-importable
dataset. This value only ever showed up in terminal/.out output before,
with no structured record -- this keeps it alongside the other
_dataset.dat files so it survives past the log file.

Usage:
    python ~/scripts/parse_homo_lumo.py --proper proper.out --output homo_lumo_dataset.dat
"""

import argparse
import re

BLOCK_RE = re.compile(
    r"HOMO\s*:\s*[-\d.]+\s*H\s*=\s*([-\d.]+)\s*eV\s*"
    r"LUMO\s*:\s*[-\d.]+\s*H\s*=\s*([-\d.]+)\s*eV\s*"
    r"HOMO-LUMO gap:\s*[-\d.]+\s*H\s*=\s*\+?([-\d.]+)\s*eV",
    re.MULTILINE,
)


def parse_homo_lumo(path):
    with open(path, "r") as f:
        text = f.read()

    match = BLOCK_RE.search(text)
    if not match:
        raise ValueError(f"No 'HOMO-LUMO Separation' block matched in {path}.")

    homo_ev, lumo_ev, gap_ev = (float(x) for x in match.groups())
    return {"homo_eV": homo_ev, "lumo_eV": lumo_ev, "gap_eV": gap_ev}


def main():
    parser = argparse.ArgumentParser(
        description="Parse Turbomole HOMO-LUMO gap into a Veusz-importable dataset.")
    parser.add_argument("--proper", required=True,
                         help="Path to a -proper run's .out file (e.g. proper.out)")
    parser.add_argument("--output", default="homo_lumo_dataset.dat")
    args = parser.parse_args()

    values = parse_homo_lumo(args.proper)
    print(f"HOMO = {values['homo_eV']:.5f} eV, LUMO = {values['lumo_eV']:.5f} eV, "
          f"gap = {values['gap_eV']:.5f} eV")

    with open(args.output, "w") as f:
        f.write("descriptor homo_eV,lumo_eV,gap_eV\n")
        f.write(f"{values['homo_eV']:.6f} {values['lumo_eV']:.6f} {values['gap_eV']:.6f}\n")

    print(f"Veusz dataset saved -> {args.output}")


if __name__ == "__main__":
    main()
