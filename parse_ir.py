#!/usr/bin/env python3
"""
parse_ir.py

Parses the Turbomole vibspectrum file (aoforce output) into a
Veusz-importable dataset, mirroring Fig. 7 of Januario & Cabral (2026).

The first 6 modes (near-zero wavenumber) correspond to translations
and rotations, not real vibrations, and are excluded from the output.
Raw IR intensities (km/mol) are kept as-is -- no normalization applied,
consistent with keeping the data as unfiltered as possible.

Usage:
    python ~/scripts/parse_ir.py --vibspectrum vibspectrum --output ir_dataset.dat
"""

import argparse
import re

ROW_RE = re.compile(
    r"^\s*(\d+)\s+(?:([a-zA-Z0-9']+)\s+)?([-\d.]+)\s+([-\d.]+)\s+(\S+)\s+(\S+)\s*$"
)


def parse_vibspectrum(path):
    modes = []
    with open(path, "r") as f:
        for line in f:
            if line.strip().startswith("$") or line.strip().startswith("#"):
                continue
            match = ROW_RE.match(line)
            if match:
                mode, symmetry, wavenumber, intensity, ir_flag, raman_flag = match.groups()
                modes.append({
                    "mode": int(mode),
                    "symmetry": symmetry if symmetry else "",
                    "wavenumber_cm1": float(wavenumber),
                    "ir_intensity_kmmol": float(intensity),
                    "ir_active": ir_flag.strip().upper() == "YES",
                    "raman_active": raman_flag.strip().upper() == "YES",
                })

    if not modes:
        raise ValueError(f"No vibrational mode rows matched in {path}. "
                          f"Check the file format.")

    return modes


def main():
    parser = argparse.ArgumentParser(
        description="Parse Turbomole vibspectrum into a Veusz-importable IR dataset.")
    parser.add_argument("--vibspectrum", required=True, help="Path to vibspectrum file")
    parser.add_argument("--output", default="ir_dataset.dat",
                         help="Output .dat filename (Veusz descriptor format)")
    parser.add_argument("--exclude-first", type=int, default=6,
                         help="Number of leading (translational/rotational) modes to exclude")
    args = parser.parse_args()

    all_modes = parse_vibspectrum(args.vibspectrum)
    print(f"Parsed {len(all_modes)} total modes from {args.vibspectrum}")

    vibrational_modes = [m for m in all_modes if m["mode"] > args.exclude_first]
    print(f"Kept {len(vibrational_modes)} vibrational modes "
          f"(excluded the first {args.exclude_first} as translation/rotation)")
    print(f"Wavenumber range: {vibrational_modes[0]['wavenumber_cm1']:.2f} to "
          f"{vibrational_modes[-1]['wavenumber_cm1']:.2f} cm-1")

    with open(args.output, "w") as f:
        f.write("descriptor mode,wavenumber_cm1,ir_intensity_kmmol\n")
        for m in vibrational_modes:
            f.write(f"{m['mode']} {m['wavenumber_cm1']:.4f} {m['ir_intensity_kmmol']:.6f}\n")

    print(f"Veusz dataset saved -> {args.output}")


if __name__ == "__main__":
    main()