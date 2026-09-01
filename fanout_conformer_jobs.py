#!/usr/bin/env python3
"""
fanout_conformer_jobs.py

Given a Turbomole r2SCAN-3c template folder (control+basis+mos, from
setup_r2scan3c_template.py) and a multi-frame ensemble .xyz (e.g. Layer
A's population_filtered.xyz), creates one job directory per frame:
control/basis/mos copied from the template, 'coord' generated from
that frame's own geometry. Each directory is ready for `jobex` --
independent of the others, meant to be driven by a SLURM job array.

Usage:
    python ~/work_turbomole/scripts/fanout_conformer_jobs.py \
        --template template_dir --ensemble population_filtered.xyz \
        --output-dir conformer_jobs
"""

import argparse
import os
import shutil
import subprocess


def read_xyz_frames(path):
    with open(path) as f:
        lines = f.readlines()
    frames = []
    i = 0
    while i < len(lines):
        natoms = int(lines[i].split()[0])
        frames.append("".join(lines[i:i + natoms + 2]))
        i += natoms + 2
    return frames


def main():
    parser = argparse.ArgumentParser(
        description="Fan out one Turbomole r2SCAN-3c template into one "
                    "job directory per conformer in an ensemble .xyz.")
    parser.add_argument("--template", required=True,
                         help="Template dir with control/basis/mos "
                              "(from setup_r2scan3c_template.py)")
    parser.add_argument("--ensemble", required=True,
                         help="Multi-frame .xyz, e.g. layer_a/population_filtered.xyz")
    parser.add_argument("--output-dir", required=True,
                         help="Parent folder to create job_0001/, job_0002/, ... in")
    args = parser.parse_args()

    for fname in ("control", "basis", "mos"):
        if not os.path.isfile(os.path.join(args.template, fname)):
            raise SystemExit(f"Template '{args.template}' is missing '{fname}'.")

    frames = read_xyz_frames(args.ensemble)
    if os.path.exists(args.output_dir):
        raise FileExistsError(f"{args.output_dir} already exists -- refusing to overwrite.")
    os.makedirs(args.output_dir)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    for i, frame in enumerate(frames, start=1):
        job_dir = os.path.join(args.output_dir, f"job_{i:04d}")
        os.makedirs(job_dir)
        for fname in ("control", "basis", "mos"):
            shutil.copy(os.path.join(args.template, fname), os.path.join(job_dir, fname))
        xyz_path = os.path.join(job_dir, "structure.xyz")
        with open(xyz_path, "w") as f:
            f.write(frame)
        subprocess.run(
            ["python3", os.path.join(script_dir, "xyz_to_coord.py"),
             "--xyz", xyz_path, "--output", os.path.join(job_dir, "coord"), "--force"],
            check=True, capture_output=True,
        )

    print(f"Created {len(frames)} job directories under {args.output_dir}/ "
          f"(job_0001 .. job_{len(frames):04d})")


if __name__ == "__main__":
    main()
