#!/usr/bin/env python3
"""
Render Artepillin C's DFT-relaxed geometry (coord.xyz) as a VESTA-style
ray-traced ball-and-stick figure via ASE + POV-Ray.

Note: the conda-forge povray 3.7.0.10 build in this env has a parser
bug -- any camera block using ASE's default right/up/direction vector
style (no explicit `angle`) throws "Viewing angle has to be smaller
than 180 degrees" regardless of the actual implied angle. Fix: replace
the camera block with an explicit perspective camera (location + angle
+ unit right/up), placed far enough away that perspective distortion
stays small (near-orthographic look).

Usage:
    ~/miniconda3/envs/vasp_env/bin/python3 plot_molecule_vesta_style.py \
        --xyz coord.xyz --labels bond_order_dataset_labels.csv \
        --output molecule_structure_vesta
"""
import argparse
import csv
import re
import subprocess

import numpy as np
from ase.io import read, write

REAL_BOND_SEN_THRESHOLD = 0.5
RADII = {"C": 0.40, "O": 0.42, "H": 0.25}


def read_real_bonds(path):
    bonds = []
    with open(path) as f:
        for row in csv.DictReader(f):
            if float(row["shared_electron_number"]) > REAL_BOND_SEN_THRESHOLD:
                bonds.append((int(row["atom1_index"]) - 1, int(row["atom2_index"]) - 1))
    return bonds


def fix_camera_block(pov_path, camera_dist, angle_deg):
    txt = open(pov_path).read()
    new_cam = (
        "camera {\n"
        "  perspective\n"
        f"  location <0, 0, {camera_dist:.2f}>\n"
        "  look_at <0, 0, 0>\n"
        f"  angle {angle_deg:.2f}\n"
        "  right <1.42, 0, 0>\n"
        "  up <0, 1, 0>\n"
        "}"
    )
    txt2, n = re.subn(r"camera \{.*?look_at <0,0,0>\}", new_cam, txt,
                       count=1, flags=re.DOTALL)
    if n != 1:
        raise RuntimeError("Could not locate camera block to patch")
    open(pov_path, "w").write(txt2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--output", default="molecule_structure_vesta")
    ap.add_argument("--canvas-width", type=int, default=1200)
    args = ap.parse_args()

    atoms = read(args.xyz)
    atoms.translate(-atoms.get_center_of_mass())
    radii = np.array([RADII[s] for s in atoms.get_chemical_symbols()])

    bondatoms = read_real_bonds(args.labels)

    pos = atoms.get_positions()
    R = np.linalg.norm(pos, axis=1).max() + radii.max()
    camera_dist = 6 * R
    angle_deg = 2 * np.degrees(np.arctan((R * 1.25) / camera_dist))

    pov_path = f"{args.output}.pov"
    ini_path = f"{args.output}.ini"
    write(pov_path, atoms, format="pov", radii=radii,
          povray_settings=dict(canvas_width=args.canvas_width,
                                bondlinewidth=0.12, bondatoms=bondatoms))

    fix_camera_block(pov_path, camera_dist, angle_deg)

    import sys, os
    povray_bin = os.path.join(os.path.dirname(sys.executable), "povray")
    subprocess.run([povray_bin, ini_path], check=True)
    print(f"Saved -> {args.output}.png")


if __name__ == "__main__":
    main()
