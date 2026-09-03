#!/usr/bin/env python3
"""Renders a static isosurface image from a Gaussian .cub file (headless, no GUI/X11).

Usage:
    python render_cube_isosurface.py --cube 81a.cub --out homo.png --isovalue 0.03 --title "HOMO (81a)"
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from skimage import measure

BOHR_TO_ANG = 0.529177210903

ELEMENT_STYLE = {
    1: dict(color="#dddddd", size=25, label="H"),
    6: dict(color="#404040", size=55, label="C"),
    7: dict(color="#3050f0", size=55, label="N"),
    8: dict(color="#e00000", size=55, label="O"),
}

BOND_CUTOFF_ANG = 1.75


def read_cube(path):
    with open(path) as fh:
        lines = fh.readlines()

    natoms, origin = int(lines[2].split()[0]), np.array(lines[2].split()[1:4], dtype=float)
    n = np.zeros(3, dtype=int)
    step = np.zeros(3)
    for i in range(3):
        parts = lines[3 + i].split()
        n[i] = int(parts[0])
        step[i] = float(parts[1 + i])  # assumes axis-aligned (orthorhombic) grid

    atom_lines = lines[6:6 + natoms]
    atoms = []
    for al in atom_lines:
        p = al.split()
        atoms.append((int(p[0]), float(p[2]), float(p[3]), float(p[4])))

    data_lines = lines[6 + natoms:]
    values = np.fromstring(" ".join(data_lines), sep=" ")
    values = values.reshape(n[0], n[1], n[2])

    return dict(origin=origin, step=step, n=n, atoms=atoms, values=values)


def add_isosurface(ax, values, spacing, origin, level, color, alpha):
    if level > 0:
        mask_present = np.any(values >= level)
    else:
        mask_present = np.any(values <= level)
    if not mask_present:
        return False

    verts, faces, _, _ = measure.marching_cubes(values, level=level, spacing=spacing)
    verts = verts + origin
    mesh = Poly3DCollection(verts[faces], alpha=alpha)
    mesh.set_facecolor(color)
    mesh.set_edgecolor("none")
    ax.add_collection3d(mesh)
    return True


def add_atoms_and_bonds(ax, atoms):
    coords = np.array([(x, y, z) for (_, x, y, z) in atoms]) * BOHR_TO_ANG
    elements = [a[0] for a in atoms]

    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            d = np.linalg.norm(coords[i] - coords[j])
            if d <= BOND_CUTOFF_ANG:
                ax.plot(*zip(coords[i], coords[j]), color="#888888", linewidth=1.2, zorder=1)

    for elem in set(elements):
        style = ELEMENT_STYLE.get(elem, dict(color="#cccccc", size=40, label=str(elem)))
        idx = [k for k, e in enumerate(elements) if e == elem]
        ax.scatter(coords[idx, 0], coords[idx, 1], coords[idx, 2],
                   c=style["color"], s=style["size"], edgecolor="black",
                   linewidth=0.4, depthshade=False, zorder=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cube", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--isovalue", type=float, default=0.03,
                     help="isosurface value in e^-1/2 bohr^-3/2; negative lobe uses the same magnitude")
    ap.add_argument("--title", default=None)
    ap.add_argument("--pos-color", default="#3060ff")
    ap.add_argument("--neg-color", default="#ff5030")
    ap.add_argument("--alpha", type=float, default=0.65)
    ap.add_argument("--elev", type=float, default=20)
    ap.add_argument("--azim", type=float, default=-60)
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()

    cube = read_cube(args.cube)
    spacing = tuple(cube["step"] * BOHR_TO_ANG)
    origin = cube["origin"] * BOHR_TO_ANG

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")

    has_pos = add_isosurface(ax, cube["values"], spacing, origin, args.isovalue, args.pos_color, args.alpha)
    has_neg = add_isosurface(ax, cube["values"], spacing, origin, -args.isovalue, args.neg_color, args.alpha)
    if not (has_pos or has_neg):
        raise SystemExit(f"No isosurface found at +/-{args.isovalue} for {args.cube} "
                          f"(data range: [{cube['values'].min():.4f}, {cube['values'].max():.4f}]) "
                          "-- lower --isovalue.")

    add_atoms_and_bonds(ax, cube["atoms"])

    coords = np.array([(x, y, z) for (_, x, y, z) in cube["atoms"]]) * BOHR_TO_ANG
    center = coords.mean(axis=0)
    span = (coords.max(axis=0) - coords.min(axis=0)).max() / 2 + 2.0
    ax.set_xlim(center[0] - span, center[0] + span)
    ax.set_ylim(center[1] - span, center[1] + span)
    ax.set_zlim(center[2] - span, center[2] + span)
    ax.set_box_aspect([1, 1, 1])
    ax.set_axis_off()
    ax.view_init(elev=args.elev, azim=args.azim)

    if args.title:
        ax.set_title(f"{args.title}\nisovalue = {args.isovalue} e$^{{-1/2}}$ bohr$^{{-3/2}}$")
    else:
        ax.set_title(f"isovalue = {args.isovalue} e$^{{-1/2}}$ bohr$^{{-3/2}}$")

    fig.tight_layout()
    fig.savefig(args.out, dpi=args.dpi, transparent=False, facecolor="white")
    print(f"wrote {args.out} (isovalue={args.isovalue}, atoms={len(cube['atoms'])}, "
          f"grid={cube['n'][0]}x{cube['n'][1]}x{cube['n'][2]})")


if __name__ == "__main__":
    main()
