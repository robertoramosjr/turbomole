#!/usr/bin/env python3
"""
setup_population_job.py

Builds a population-analysis job directory (bond order/SEN via
`dscf -proper`, plus a density cube for Bader via `ridft -proper`) from
an already-converged HSE06/def2-TZVP conformer directory, reusing its
converged `mos` (no new SCF from scratch, no new Hessian).

Ensures `$pop mulliken loewdin nbo paboon` is present in the copied
`control` (population analysis, for parse_bond_order.py /
parse_mulliken.py -- often already there from `define`, added only if
missing), and writes the per-conformer `$pointval fmt=cub` block
(3D real-space grid for the Bader cube) to a SEPARATE file
(`pointval_block.txt`) rather than adding it to `control` directly.

Why separate: confirmed empirically that with `$pop` and `$pointval`
both present in the same `control`, EITHER `dscf -proper` or `ridft
-proper` will attempt the expensive (~100s) grid evaluation, since
`-proper` honors `$pointval` regardless of which binary reads it --
running both steps against one control this way silently duplicates
the cube computation instead of splitting cheap (pop, ~20s) from
expensive (cube, ~100s) work. Paper 1 avoided this by keeping its
cube-generating run in a separate `badergrid/` control from the plain
`dscf -proper` one; the job script here instead runs `dscf -proper`
first (control without `$pointval`), THEN appends `pointval_block.txt`
to `control` and runs `ridft -proper` for the cube.

Grid box: conformers here have different atomic coordinates
(CREST-sampled geometries, unlike Paper 1's single reference), so the
box is computed PER CONFORMER from its own coord file: 5.0 bohr
padding beyond the atomic bounding box on every side, ~0.15 bohr point
spacing (both values reverse-engineered from Paper 1's own grid: its
box was exactly 5.0 bohr beyond the atomic extent in x/y/z, and
points-per-span works out to ~0.15 bohr spacing in all three
dimensions).

Usage:
    python ~/work_turbomole/scripts/setup_population_job.py \
        --source-dir layer_b/labeled_jobs/job_0011 \
        --output-dir layer_b/population_jobs/job_0011
"""

import argparse
import os
import shutil

PADDING_BOHR = 5.0
GRID_SPACING_BOHR = 0.15


def read_coord_bounds(coord_path):
    xs, ys, zs = [], [], []
    with open(coord_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 4:
                try:
                    x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                except ValueError:
                    continue
                xs.append(x)
                ys.append(y)
                zs.append(z)
    return (min(xs), max(xs)), (min(ys), max(ys)), (min(zs), max(zs))


def pointval_block(coord_path):
    (xmin, xmax), (ymin, ymax), (zmin, zmax) = read_coord_bounds(coord_path)
    bounds = []
    for lo, hi in [(xmin, xmax), (ymin, ymax), (zmin, zmax)]:
        lo_p, hi_p = lo - PADDING_BOHR, hi + PADDING_BOHR
        npoints = round((hi_p - lo_p) / GRID_SPACING_BOHR)
        bounds.append((lo_p, hi_p, npoints))
    (gx0, gx1, nx), (gy0, gy1, ny), (gz0, gz1, nz) = bounds
    return (
        "$pointval fmt=cub\n"
        f"   grid1 vector 1 0 0 range {gx0:.3f},{gx1:.3f} points {nx}\n"
        f"   grid2 vector 0 1 0 range {gy0:.3f},{gy1:.3f} points {ny}\n"
        f"   grid3 vector 0 0 1 range {gz0:.3f},{gz1:.3f} points {nz}\n"
        "   origin 0 0 0\n"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Build a population-analysis (bond order + Bader cube) "
                    "job directory from an already-converged conformer.")
    parser.add_argument("--source-dir", required=True,
                         help="Already-converged HSE06 job dir (control+basis+"
                              "auxbasis+mos+coord)")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    if os.path.exists(args.output_dir):
        raise FileExistsError(f"{args.output_dir} already exists -- refusing to overwrite.")
    os.makedirs(args.output_dir)

    for fname in ("control", "basis", "auxbasis", "mos", "coord"):
        src = os.path.join(args.source_dir, fname)
        if os.path.isfile(src):
            shutil.copy(src, os.path.join(args.output_dir, fname))

    control_path = os.path.join(args.output_dir, "control")
    with open(control_path) as f:
        text = f.read()

    # $pointval is intentionally NOT added here. Confirmed empirically:
    # with $pop and $pointval both present in the same control, EITHER
    # `dscf -proper` or `ridft -proper` will attempt the (expensive,
    # ~100s) grid/cube evaluation, since `-proper` honors whatever
    # $pointval says regardless of which binary reads it -- running
    # both steps against one control this way silently duplicates the
    # cube computation instead of splitting cheap (pop) from expensive
    # (cube) work. Paper 1 avoided this by keeping the cube-generating
    # run in a separate `badergrid/` control from the plain `dscf
    # -proper` one; this project's job script instead injects
    # $pointval into a copy of control just before the `ridft -proper`
    # step, after `dscf -proper` (bond order/SEN, no $pointval, fast)
    # has already run and written its own pop.out.
    addition = ""
    if "$pop " not in text and "$pop\n" not in text:
        addition += "$pop mulliken loewdin nbo paboon\n"
    if addition:
        text = text.replace("$end", addition + "$end", 1)
        with open(control_path, "w") as f:
            f.write(text)

    grid_path = os.path.join(args.output_dir, "pointval_block.txt")
    with open(grid_path, "w") as f:
        f.write(pointval_block(os.path.join(args.output_dir, "coord")))

    print(f"Population job ready in {args.output_dir}/ "
          f"(pointval_block.txt has the per-conformer cube grid, to be "
          f"injected into control by the job script AFTER dscf -proper runs)")


if __name__ == "__main__":
    main()
