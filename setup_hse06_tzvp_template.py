#!/usr/bin/env python3
"""
setup_hse06_tzvp_template.py

Builds a production-level HSE06-D3(BJ)/def2-TZVP Turbomole template
(control + basis + auxbasis + mos) for one protonation state -- the
"expensive target" level for Layer B's B1 (aoforce/IR) and B2
(escf/TD-DFT) surrogates. Same functional/basis/dispersion/grid as
Paper 1's production conformer, confirmed against its `control`
(functional hse06, gridsize m4, $senex, $scfconv 8, $rij, $disp3 bj) --
this is what keeps the two papers comparable, not a fresh choice.

Unlike PBE0 (Layer A's IR/UV-Vis functional), HSE06 is a *screened*-
exchange hybrid and only needs RI-J (not RI-K) to be tractable --
confirmed already in prepare_functional_benchmark.py's docstring ("a
plain-RI-J PBE0 single point did not finish one SCF cycle in 3
minutes... HSE06 -proper finished in 35s"). def2-TZVP is a mainstream,
fully-covered basis (unlike r2SCAN-3c's def2-mTZVP), so no auxiliary-
basis surprises expected here.

This builds a template for property calculations (aoforce Hessian +
TD-DFT) on top of an already-optimized geometry (Layer A's r2SCAN-3c
structures) -- it does NOT reoptimize geometry at this level. Full
production reoptimization (matching Paper 1's tight $statpt thresholds)
is reserved for Estagio 6's final top-N selection, not this labeling
step.

Usage:
    python ~/work_turbomole/scripts/setup_hse06_tzvp_template.py \
        --xyz structure.xyz --charge -1 --nstates 60 \
        --output-dir template_dir
"""

import argparse
import os
import subprocess
import sys
import time

EXTRA_BLOCK = (
    "$senex\n"
    "$rij\n"
    "$disp3 bj\n"
    "$scfconv 8\n"
    "$denconv 1d-7\n"
    "$scfinstab rpas\n"
    "$soes\n"
    " a           {nstates}\n"
    "$rpacor 500\n"
)


def advance_to(child, patterns, max_steps=25, step_wait=3):
    import pexpect
    for _ in range(max_steps):
        idx = child.expect(patterns + [pexpect.TIMEOUT], timeout=step_wait)
        if idx < len(patterns):
            return idx
        child.sendline("")
    raise RuntimeError(f"define: did not reach any of {patterns}")


def run_define(work_dir, charge, log_path):
    import pexpect
    child = pexpect.spawn("define", cwd=work_dir, timeout=20, encoding="utf-8")
    with open(log_path, "w") as logf:
        child.logfile = logf

        def send(s, wait=0.8):
            child.sendline(s)
            time.sleep(wait)

        send("")
        send("")
        send("a coord")
        send("*")
        send("no")
        send("b all def2-TZVP", wait=2.0)
        send("*", wait=1.0)
        send("eht")
        send("")
        send(str(charge))
        send("")
        send("dft")
        send("func hse06")
        send("grid m4")
        send("on")
        send("")
        send("*", wait=2.0)
        advance_to(child, ["GENERAL MENU"], max_steps=5, step_wait=3)
        send("*", wait=2.0)
        try:
            child.expect(pexpect.EOF, timeout=25)
        finally:
            child.close()


def main():
    parser = argparse.ArgumentParser(
        description="Build a reusable HSE06-D3(BJ)/def2-TZVP Turbomole "
                    "template (production level) for one protonation state.")
    parser.add_argument("--xyz", required=True)
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--nstates", type=int, default=60,
                         help="TD-DFT singlet roots [default: 60, matches Paper 1]")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    try:
        import pexpect  # noqa: F401
    except ImportError:
        sys.exit("This script requires pexpect: pip install pexpect")

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
    run_define(args.output_dir, args.charge, log_path)

    control_path = os.path.join(args.output_dir, "control")
    if not os.path.isfile(control_path):
        raise RuntimeError(f"define did not produce a control file -- check {log_path}")

    with open(control_path) as f:
        text = f.read()
    if "hse06" not in text:
        raise RuntimeError(
            f"control does not show the expected hse06 setup -- "
            f"check {log_path} and {control_path} by hand.")

    # define's own default $scfconv would otherwise duplicate (and
    # potentially conflict with) the $scfconv 8 added below to match
    # Paper 1 exactly.
    import re
    text = re.sub(r"^\$scfconv\s+\d+\s*\n", "", text, count=1, flags=re.MULTILINE)

    text = text.replace("$end", EXTRA_BLOCK.format(nstates=args.nstates) + "$end", 1)
    with open(control_path, "w") as f:
        f.write(text)

    print(f"Template ready in {args.output_dir}/ (control, basis, auxbasis, mos, "
          f"coord; charge={args.charge}, HSE06-D3(BJ)/def2-TZVP, grid m4, "
          f"RI-J, senex, TD-DFT {args.nstates} singlet states -- Paper 1 "
          f"production settings).")


if __name__ == "__main__":
    main()
