#!/usr/bin/env python3
"""
poscar_to_xyz.py

Converts a VASP POSCAR (Direct/fractional coordinates, orthogonal cell
only) into a standard XYZ file in Angstrom, for downstream conversion
to a Turbomole 'coord' file via xyz_to_coord.py. Used to compare
VASP-side relaxations against Turbomole production runs on the same
starting geometry.

Only supports orthogonal (diagonal) lattice vectors and 'Direct'
(fractional) coordinate mode -- VASP's most common isolated-molecule-
in-a-box convention. Cartesian-mode POSCARs or non-orthogonal cells
are not handled.

Usage:
    python poscar_to_xyz.py --poscar POSCAR --output molecule.xyz
"""

import argparse


def read_poscar(poscar_path):
    with open(poscar_path) as f:
        lines = [line.rstrip("\n") for line in f]

    comment = lines[0]
    scale = float(lines[1].split()[0])
    lattice = []
    for i in range(2, 5):
        vec = [float(x) for x in lines[i].split()]
        lattice.append(vec)
    for i, vec in enumerate(lattice):
        if abs(vec[i]) < 1e-9 or any(abs(vec[j]) > 1e-9 for j in range(3) if j != i):
            raise ValueError(
                f"poscar_to_xyz.py only supports orthogonal (diagonal) "
                f"lattice vectors -- lattice line {i} is {vec}"
            )
    box = [lattice[i][i] * scale for i in range(3)]

    elements = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    mode_line = lines[7].strip().lower()
    if not mode_line.startswith("d"):
        raise ValueError(
            f"poscar_to_xyz.py only supports 'Direct' coordinate mode, "
            f"got: {lines[7]!r}"
        )

    atoms = []
    line_idx = 8
    for element, count in zip(elements, counts):
        for _ in range(count):
            parts = lines[line_idx].split()
            fx, fy, fz = float(parts[0]), float(parts[1]), float(parts[2])
            x = fx * box[0]
            y = fy * box[1]
            z = fz * box[2]
            atoms.append((element, x, y, z))
            line_idx += 1

    return comment, atoms


def write_xyz(comment, atoms, output_path):
    with open(output_path, "w") as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{comment}\n")
        for symbol, x, y, z in atoms:
            f.write(f"{symbol:2s} {x:20.14f} {y:20.14f} {z:20.14f}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert a VASP POSCAR (Direct, orthogonal cell) to XYZ (Angstrom)."
    )
    parser.add_argument("--poscar", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    comment, atoms = read_poscar(args.poscar)
    write_xyz(comment, atoms, args.output)
    print(f"Wrote {len(atoms)} atoms to {args.output}")


if __name__ == "__main__":
    main()
