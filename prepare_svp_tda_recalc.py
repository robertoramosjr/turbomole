#!/usr/bin/env python3
"""
prepare_svp_tda_recalc.py

Selects a random sample of Wigner-sampled points from a Newton-X
`initial_condition.<nis>.<nfs>` file and builds one single-point folder
per selected point, reusing a validated Turbomole control template
(basis/method/RI/dispersion already fixed -- only the geometry changes
per folder). Intended for recalibrating a laser-window definition at a
different level of theory than the one used for the original ensemble
(e.g. def2-TZVP/RPA -> def2-SVP/TDA), on a representative subsample
rather than the full 1000 points.

The random sample uses a fixed, recorded seed so the same 100 points
can be re-identified later against the original dataset (e.g.
wigner_1000_states_dump.csv) for a point-by-point energy-shift
comparison.

Usage:
    python ~/scripts/prepare_svp_tda_recalc.py \
        --initial-condition initial_condition.1.2 \
        --n-points 100 --total-points 1000 --seed 42 \
        --template-dir svp_tda_template \
        --output-dir svp_tda_recalc
"""

import argparse
import os
import random
import shutil

BOHR_LABEL_RE_COLUMNS = 5  # label, Z, x, y, z, mass in the NX geom block


def parse_initial_condition_geometries(path):
    """Returns a list of geometry blocks (list of (label, x, y, z) tuples
    in bohr), in file order: index 0 is the equilibrium geometry, indices
    1..N are the Wigner-sampled points in generation order."""
    geometries = []
    current = None
    with open(path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("Geometry in COLUMBUS"):
                current = []
                continue
            if stripped.startswith("Velocity in NX input format"):
                if current is not None:
                    geometries.append(current)
                current = None
                continue
            if current is not None and stripped:
                parts = stripped.split()
                label, _z, x, y, z, _mass = parts[:6]
                current.append((label, float(x), float(y), float(z)))
    if not geometries:
        raise ValueError(f"No geometry blocks found in {path}")
    return geometries


def write_turbomole_coord(path, atoms):
    with open(path, "w") as f:
        f.write("$coord\n")
        for label, x, y, z in atoms:
            f.write(f"{x:20.14f}{y:20.14f}{z:20.14f}  {label}\n")
        f.write("$end\n")


def main():
    parser = argparse.ArgumentParser(
        description="Select a random sample of Wigner points and build one "
                    "single-point Turbomole folder per point from a template.")
    parser.add_argument("--initial-condition", required=True,
                         help="Path to a Newton-X initial_condition.<nis>.<nfs> file")
    parser.add_argument("--n-points", type=int, required=True,
                         help="Number of points to sample")
    parser.add_argument("--total-points", type=int, required=True,
                         help="Total number of Wigner-sampled points available "
                              "(excludes the equilibrium geometry block)")
    parser.add_argument("--seed", type=int, required=True,
                         help="Random seed for the sample (recorded in the manifest)")
    parser.add_argument("--template-dir", required=True,
                         help="Folder with a validated control/basis/auxbasis/mos "
                              "(geometry-independent; --coord in control is overwritten)")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--manifest", default="sample_manifest.csv",
                         help="Output CSV recording which point_id went into which folder")
    args = parser.parse_args()

    if os.path.exists(args.output_dir):
        raise FileExistsError(f"{args.output_dir} already exists -- refusing to overwrite.")

    geometries = parse_initial_condition_geometries(args.initial_condition)
    n_available = len(geometries) - 1  # exclude the equilibrium block (index 0)
    if n_available != args.total_points:
        raise ValueError(
            f"{args.initial_condition} has {n_available} Wigner points, "
            f"expected --total-points {args.total_points}")

    rng = random.Random(args.seed)
    selected_point_ids = sorted(rng.sample(range(1, args.total_points + 1), args.n_points))

    os.makedirs(args.output_dir)
    template_files = [f for f in os.listdir(args.template_dir)
                       if f not in ("coord",)]

    manifest_rows = []
    for point_id in selected_point_ids:
        folder = os.path.join(args.output_dir, f"P{point_id}")
        os.makedirs(folder)
        for fname in template_files:
            shutil.copy(os.path.join(args.template_dir, fname), os.path.join(folder, fname))
        write_turbomole_coord(os.path.join(folder, "coord"), geometries[point_id])
        manifest_rows.append({"point_id": point_id, "folder": f"P{point_id}"})

    with open(args.manifest, "w") as f:
        f.write("point_id,folder,seed,n_points,total_points,source\n")
        for r in manifest_rows:
            f.write(f"{r['point_id']},{r['folder']},{args.seed},{args.n_points},"
                     f"{args.total_points},{args.initial_condition}\n")

    print(f"Sampled {args.n_points} of {args.total_points} points (seed={args.seed})")
    print(f"Point IDs: {selected_point_ids}")
    print(f"Folders created under -> {args.output_dir}")
    print(f"Manifest saved -> {args.manifest}")


if __name__ == "__main__":
    main()
