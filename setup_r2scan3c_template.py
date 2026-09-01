#!/usr/bin/env python3
"""
setup_r2scan3c_template.py

Drives Turbomole's `define` (via pexpect, over a real PTY -- define's
prompts behave differently under a plain pipe) to build one r2SCAN-3c
template (control + basis + mos) for a given protonation state. Meant
to be run ONCE per state (charge), not once per conformer: the basis
assignment only depends on atom composition/charge, not geometry, so
every population-filtered conformer of that state can reuse this one
template and just swap in its own 'coord'.

r2SCAN-3c is a Grimme-group composite method: functional r2SCAN + a
def2-mTZVP basis + D4 dispersion + a matching gCP basis-set-
superposition correction. Turbomole's ridft/jobex auto-detect the
'r2scan-3c' $dft functional keyword and apply D4 + gCP automatically
at runtime, but do NOT auto-correct the basis set -- def2-mTZVP has to
be assigned explicitly during `define` (confirmed empirically: without
it, define silently keeps whatever generic default basis -- def-SV(P)
here -- was auto-assigned before the dft menu is touched, which is the
wrong, much smaller basis for this method).

Usage:
    python ~/work_turbomole/scripts/setup_r2scan3c_template.py \
        --xyz structure.xyz --charge -1 --output-dir template_dir
"""

import argparse
import os
import subprocess
import sys
import time

DEFINE_STEPS = [
    "",                      # hit return: skip reading defaults from another control
    "",                      # title: accept default
    "a coord",               # add geometry from the local 'coord' file
    "*",                     # terminate geometry specification
    "no",                    # do not use internal coordinates (cartesian, matches Paper 1)
    "b all def2-mTZVP",      # r2SCAN-3c's prescribed basis (not auto-assigned by define)
    "*",                     # accept basis, leave atomic attribute menu
    "eht",                   # extended Hueckel initial guess
    "",                      # accept default EHT parameters
    None,                    # placeholder for charge, filled in at call time
    "",                      # accept default MO occupation
    "dft",                   # enter DFT menu
    "func r2scan-3c",        # composite method (drives auto D4 + gCP at runtime)
    "grid m4",               # match Paper 1's production grid choice
    "on",                    # switch DFT on
    "",                      # leave DFT menu
    "*",                     # end define session, write control
]


def run_define(work_dir, charge, log_path):
    import pexpect
    steps = [charge if s is None else s for s in DEFINE_STEPS]
    child = pexpect.spawn("define", cwd=work_dir, timeout=20, encoding="utf-8")
    with open(log_path, "w") as logf:
        child.logfile = logf
        for i, step in enumerate(steps):
            wait = 2.5 if "b all" in step else 0.8
            child.sendline(step)
            time.sleep(wait)
        try:
            child.expect(pexpect.EOF, timeout=20)
        finally:
            child.close()


def main():
    parser = argparse.ArgumentParser(
        description="Build a reusable r2SCAN-3c/def2-mTZVP Turbomole "
                    "template (control+basis+mos) for one protonation "
                    "state via a scripted `define` session.")
    parser.add_argument("--xyz", required=True, help="Representative .xyz for this state")
    parser.add_argument("--charge", type=int, required=True,
                         help="Total molecular charge for this protonation state")
    parser.add_argument("--output-dir", required=True,
                         help="Folder to create with the finished template")
    args = parser.parse_args()

    try:
        import pexpect  # noqa: F401
    except ImportError:
        sys.exit("This script requires pexpect (define's prompts need a real "
                  "PTY): pip install pexpect")

    if os.path.exists(args.output_dir):
        raise FileExistsError(f"{args.output_dir} already exists -- refusing to overwrite.")
    os.makedirs(args.output_dir)

    coord_path = os.path.join(args.output_dir, "coord")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    subprocess.run(
        ["python3", os.path.join(script_dir, "xyz_to_coord.py"),
         "--xyz", os.path.abspath(args.xyz), "--output", coord_path],
        check=True,
    )

    log_path = os.path.join(args.output_dir, "define_session.log")
    run_define(args.output_dir, str(args.charge), log_path)

    control_path = os.path.join(args.output_dir, "control")
    if not os.path.isfile(control_path):
        raise RuntimeError(
            f"define did not produce a control file -- check {log_path}")
    with open(control_path) as f:
        control_text = f.read()
    if "r2scan-3c" not in control_text or "def2-mTZVP" not in open(
            os.path.join(args.output_dir, "basis")).read():
        raise RuntimeError(
            f"control/basis do not show the expected r2scan-3c/def2-mTZVP "
            f"setup -- check {log_path} and {control_path} by hand.")

    # Enable RI-J (define's default session here doesn't turn it on, and
    # ridft refuses to run without it) and add loosened Layer-A geometry
    # convergence thresholds (not Paper 1's production thresholds --
    # see filter_boltzmann_population.py / job array script docstrings
    # for why a screening stage does not need production-tight geometry
    # convergence).
    with open(control_path) as f:
        text = f.read()
    text = text.replace(
        "$end",
        "$statpt\n"
        "   threchange   1.0d-6\n"
        "   thrmaxdispl  1.0d-3\n"
        "   thrmaxgrad   1.0d-3\n"
        "   thrrmsdispl  5.0d-4\n"
        "   thrrmsgrad   5.0d-4\n"
        "$rij\n"
        "$end",
        1,
    )
    with open(control_path, "w") as f:
        f.write(text)

    print(f"Template ready in {args.output_dir}/ "
          f"(control, basis, mos, coord; charge={args.charge}, "
          f"r2scan-3c/def2-mTZVP, grid m4, RI-J on, "
          f"loosened Layer-A $statpt thresholds).")


if __name__ == "__main__":
    main()
