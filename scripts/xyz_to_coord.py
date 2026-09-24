#!/usr/bin/env python3
"""
xyz_to_coord.py

Converts a standard XYZ geometry file into a Turbomole 'coord' file
(atomic units / Bohr, lowercase element symbols, wrapped in a
'$coord' / '$end' data group). This is the first step of the pipeline,
run before 'define' / Stage 1 (ridft, escf, aoforce, ...).

Usage:
    python ~/scripts/xyz_to_coord.py --xyz molecule.xyz --output coord
"""

import argparse

BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM


def read_xyz(xyz_path):
    """Read a single-frame XYZ file, return list of (symbol, x, y, z) in Angstrom."""
    with open(xyz_path, "r") as f:
        lines = [line.rstrip("\n") for line in f]

    if len(lines) < 2:
        raise ValueError(f"'{xyz_path}' is too short to be a valid XYZ file.")

    natoms = int(lines[0].split()[0])
    atom_lines = lines[2:2 + natoms]

    if len(atom_lines) != natoms:
        raise ValueError(
            f"'{xyz_path}' declares {natoms} atoms but only "
            f"{len(atom_lines)} coordinate lines were found."
        )

    atoms = []
    for line in atom_lines:
        parts = line.split()
        symbol, x, y, z = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
        atoms.append((symbol.lower(), x, y, z))

    return atoms


def write_turbomole_coord(atoms, output_path):
    """Write atoms (symbol, x, y, z in Angstrom) as a Turbomole coord file (Bohr)."""
    with open(output_path, "w") as f:
        f.write("$coord\n")
        for symbol, x, y, z in atoms:
            xb = x * ANGSTROM_TO_BOHR
            yb = y * ANGSTROM_TO_BOHR
            zb = z * ANGSTROM_TO_BOHR
            f.write(f"{xb:20.14f}{yb:20.14f}{zb:20.14f}  {symbol}\n")
        f.write("$end\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert an XYZ geometry file into a Turbomole 'coord' file "
                    "(Angstrom -> Bohr, lowercase symbols, $coord/$end block)."
    )
    parser.add_argument("--xyz", required=True, help="Path to the input .xyz file")
    parser.add_argument("--output", default="coord",
                         help="Output filename (default: 'coord', as required by Turbomole)")
    parser.add_argument("--force", action="store_true",
                         help="Overwrite the output file if it already exists")
    args = parser.parse_args()

    import os
    if not os.path.isfile(args.xyz):
        raise SystemExit(f"Input file '{args.xyz}' not found.")

    if os.path.exists(args.output) and not args.force:
        raise SystemExit(
            f"'{args.output}' already exists. Use --force to overwrite."
        )

    atoms = read_xyz(args.xyz)
    write_turbomole_coord(atoms, args.output)

    print(f"Converted '{args.xyz}' -> '{args.output}' "
          f"({len(atoms)} atoms, Turbomole format, coordinates in Bohr).")


if __name__ == "__main__":
    main()