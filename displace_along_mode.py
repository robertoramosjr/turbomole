#!/usr/bin/env python3
"""Displace a TURBOMOLE coord file along one mode of vib_normal_modes.

Used to leave a saddle point (imaginary frequency) before re-optimizing.
The mode is rescaled so the largest single-atom displacement equals
--max-disp (Bohr); --sign flips the direction.

Usage (from the calculation folder):
    python displace_along_mode.py --coord coord --modes vib_normal_modes \
        --mode 1 --max-disp 0.4 --output coord_displaced
"""
import argparse
import math
import sys


def read_modes(path, n_cart):
    """Return {mode_index: [3N floats]} from a $vibrational normal modes file."""
    modes = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if line.startswith("$") or len(parts) < 3:
                continue
            block, row = int(parts[0]), int(parts[1])
            for k, val in enumerate(parts[2:]):
                modes.setdefault((block - 1) * 5 + k + 1, [0.0] * n_cart)[row - 1] = float(val.replace("D", "E"))
    return modes


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--coord", default="coord")
    ap.add_argument("--modes", default="vib_normal_modes")
    ap.add_argument("--mode", type=int, default=1, help="mode number as in vibspectrum")
    ap.add_argument("--max-disp", type=float, default=0.4, help="largest atomic displacement (Bohr)")
    ap.add_argument("--sign", type=int, choices=(-1, 1), default=1)
    ap.add_argument("--output", default="coord_displaced")
    args = ap.parse_args()

    lines = open(args.coord).read().splitlines()
    atom_idx = [i for i, l in enumerate(lines) if len(l.split()) == 4 and l.split()[3].isalpha()
                and not l.startswith("$")]
    n_cart = 3 * len(atom_idx)
    modes = read_modes(args.modes, n_cart)
    if args.mode not in modes:
        sys.exit(f"mode {args.mode} not found in {args.modes}")
    vec = modes[args.mode]
    peak = max(math.sqrt(sum(vec[3 * a + c] ** 2 for c in range(3))) for a in range(len(atom_idx)))
    scale = args.sign * args.max_disp / peak
    for a, i in enumerate(atom_idx):
        x, y, z, sym = lines[i].split()
        new = [float(v) + scale * vec[3 * a + c] for c, v in enumerate((x, y, z))]
        lines[i] = "%22.14f %22.14f %22.14f  %s" % (*new, sym)
    with open(args.output, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"{len(atom_idx)} atoms; mode {args.mode} scaled by {scale:.4f}; wrote {args.output}")


if __name__ == "__main__":
    main()
