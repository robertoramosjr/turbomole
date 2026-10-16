#!/usr/bin/env python3
"""
setup_cosmo_job.py

Builds a COSMO(water) single-point job directory from an
already-converged HSE06/def2-TZVP conformer directory (frozen
geometry, no reoptimization) -- for the pKa thermodynamic cycle
(protonation_and_conformer_investigation.md Sec. 9.32-9.33): G(aq) =
G(gas) + G_solv.

Drives `cosmoprep` (interactive, like `define`) automatically. Prompt
sequence mapped by direct exploration in a scratch copy (2026-09-14):
  1. epsilon -> 78.4 (water; the only prompt NOT left at its default --
     everything else below is accepted as-is)
  2. refind -> default (none)
  3. LR terms -> default (off)
  4. COSMO RF equil. -> default (not set)
  5. nppa -> default (1082)
  6. nspa -> default (92)
  7. disex -> default (10.00)
  8. rsolv -> default (min(rad(h)))
  9. routf -> default (0.85)
  10. cavity -> default (closed)
  11. amat -> default ("*")
  12. radius definition menu -> "r all b" (assign the default/optimized
      radius to every atom), then "*" to finish
  13. COSMO output file -> default (out.ccf)
  14. correlated calculation (rimp2/ricc2) -> default (n)

KNOWN BUG worked around here: cosmoprep writes "amat file=*" into
$cosmo when the amat prompt's default is accepted with "*" -- `ridft`
rejects this ("COSMO ERROR: unknown keyword in $cosmo section detected:
amat file=*"). This script deletes that line after cosmoprep runs (the
internal default applies fine without it).

Usage:
    python ~/work_turbomole/scripts/setup_cosmo_job.py \
        --source-dir layer_b/labeled_jobs/job_0011 \
        --output-dir layer_b/cosmo_jobs/job_0011 \
        --epsilon 78.4
"""

import argparse
import os
import shutil

import pexpect

FILES_TO_COPY = ("control", "basis", "auxbasis", "mos", "coord")


def run_cosmoprep(work_dir, epsilon, timeout=10):
    child = pexpect.spawn("cosmoprep", cwd=work_dir, timeout=timeout, encoding="utf-8")
    inputs = iter([
        str(epsilon), "", "", "", "", "", "", "", "", "",  # steps 1-10
        "*",            # 11: amat default
        "r all b",      # 12: assign radius to all atoms
        "*",            # 12b: finish radius definition
        "",             # 13: cosmo output file default
        "",             # 14: correlated calc default (n)
    ])
    while True:
        idx = child.expect([pexpect.EOF, pexpect.TIMEOUT])
        if idx == 0:
            break
        try:
            child.sendline(next(inputs))
        except StopIteration:
            child.sendline("")
    return "cosmoprep ended normally" in (child.before or "")


def main():
    parser = argparse.ArgumentParser(
        description="Build a COSMO(water) single-point job directory from "
                    "an already-converged HSE06/def2-TZVP conformer, no "
                    "reoptimization.")
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--epsilon", type=float, default=78.4,
                         help="Solvent dielectric constant [default: 78.4, water]")
    args = parser.parse_args()

    if os.path.exists(args.output_dir):
        raise FileExistsError(f"{args.output_dir} already exists -- refusing to overwrite.")
    os.makedirs(args.output_dir)

    for fname in FILES_TO_COPY:
        src = os.path.join(args.source_dir, fname)
        if os.path.isfile(src):
            shutil.copy(src, args.output_dir)

    control_path = os.path.join(args.output_dir, "control")
    with open(control_path) as f:
        lines = [l for l in f if not l.strip().startswith("$last step")
                 and not l.strip().startswith("$actual step")]
    with open(control_path, "w") as f:
        f.writelines(lines)

    ok = run_cosmoprep(args.output_dir, args.epsilon)
    if not ok:
        raise SystemExit(f"cosmoprep did not end normally for {args.output_dir}")

    # work around the amat file=* bug (see docstring)
    with open(control_path) as f:
        lines = [l for l in f if "amat file=*" not in l]
    with open(control_path, "w") as f:
        f.writelines(lines)

    print(f"Built COSMO job dir: {args.output_dir}")


if __name__ == "__main__":
    main()
