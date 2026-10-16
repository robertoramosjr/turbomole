#!/usr/bin/env python3
"""
setup_hse06_opt_template.py

Builds a HSE06-D3(BJ)/def2-TZVP Turbomole template ready for `jobex`
geometry optimization (ground state, no TD-DFT) with Paper 1's TIGHT
production $statpt thresholds -- for small pKa reference compounds
(protonation_and_conformer_investigation.md Sec. 9.33/9.34), whose
optimization is cheap and must be done at the SAME level as the
artepilin C single points to keep the calibration offset meaningful.

Wraps setup_hse06_tzvp_template.py (same define automation, same
functional/basis/grid/dispersion as Paper 1), then strips its escf-
specific block ($scfinstab/$soes/$rpacor -- useless for a ground-state
optimization) and appends the tight $statpt block (Paper 1: threchange
1e-7, thrmaxgrad 5e-4, thrrmsgrad 2e-4, thrmaxdispl 5e-4, thrrmsdispl
2e-4). The cycle cap is NOT set here -- it is the `-c` argument of the
`jobex` call in the job script (user directive: 800 for optimization
stages, to guarantee convergence).

Usage:
    python ~/work_turbomole/scripts/setup_hse06_opt_template.py \
        --xyz phenol.xyz --charge 0 --output-dir phenol/acid
"""

import argparse
import os
import subprocess
import sys

STATPT_BLOCK = (
    "$statpt\n"
    "   threchange 1.0d-7\n"
    "   thrmaxgrad 5.0d-4\n"
    "   thrrmsgrad 2.0d-4\n"
    "   thrmaxdispl 5.0d-4\n"
    "   thrrmsdispl 2.0d-4\n"
)
ESCF_LINES = ("$scfinstab", "$soes", "$rpacor")


def main():
    parser = argparse.ArgumentParser(
        description="Build a HSE06/def2-TZVP jobex-ready template with "
                    "Paper 1's tight $statpt thresholds.")
    parser.add_argument("--xyz", required=True)
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    # the repo was reorganized 2026-09-24 (tracked scripts moved to
    # scripts/scripts/) -- look in both places
    template_script = next(
        (os.path.join(d, "setup_hse06_tzvp_template.py")
         for d in (script_dir, os.path.join(script_dir, "scripts"))
         if os.path.isfile(os.path.join(d, "setup_hse06_tzvp_template.py"))),
        None)
    if template_script is None:
        raise FileNotFoundError("setup_hse06_tzvp_template.py not found next to "
                                "this script nor in its scripts/ subfolder")
    subprocess.run(
        [sys.executable, template_script,
         "--xyz", args.xyz, "--charge", str(args.charge), "--nstates", "1",
         "--output-dir", args.output_dir],
        check=True,
    )

    control_path = os.path.join(args.output_dir, "control")
    with open(control_path) as f:
        lines = f.readlines()
    kept = []
    skip_next = False
    for l in lines:
        if skip_next:
            skip_next = False
            continue
        if l.startswith(ESCF_LINES):
            skip_next = l.startswith("$soes")  # its " a  N" root-count line
            continue
        kept.append(l)
    text = "".join(kept)
    if "$statpt" in text:
        raise RuntimeError("control already has a $statpt block -- refusing to duplicate.")
    text = text.replace("$end", STATPT_BLOCK + "$end", 1)
    with open(control_path, "w") as f:
        f.write(text)

    print(f"jobex-ready HSE06-D3(BJ)/def2-TZVP template in {args.output_dir}/ "
          f"(charge={args.charge}, tight Paper 1 $statpt; set the cycle cap "
          f"with `jobex -c <N>`).")


if __name__ == "__main__":
    main()
