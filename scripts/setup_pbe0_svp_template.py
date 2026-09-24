#!/usr/bin/env python3
"""
setup_pbe0_svp_template.py

Drives Turbomole's `define` (via pexpect) to build one PBE0/def2-SVP
RI-K + TD-DFT template (control + basis + auxbasis + mos) for a given
protonation state -- the IR (aoforce) + UV-Vis (escf) level used
downstream of Layer A's r2SCAN-3c geometry optimization.

Why not r2SCAN-3c for IR/UV-Vis too: confirmed empirically that both
`aoforce` and `escf` abort on this Turbomole build for r2scan-3c with
"Invalid value of nfun in <mgga_r2>!" -- meta-GGA analytic
Hessian/TD-DFT response is not implemented here (the DFT-functionals.xml
FREQ=0/UVVIS=0 flags for r2-SCAN-3c were correct, not just a `define`
menu restriction as an earlier, too-hasty aoforce test suggested).

Why def2-SVP and not def2-mTZVP (r2SCAN-3c's own basis): confirmed
empirically that RI-K's auxiliary basis assignment silently produces an
empty/unusable auxbasis file for def2-mTZVP ("Problem reading basis
set(s)" from ridft) -- def2-mTZVP is a bespoke basis for the "-3c"
composite recipe, not a mainstream basis with full RI-JK library
coverage. def2-SVP is one of Turbomole's most standard, fully-covered
basis sets and works cleanly with RI-K.

PBE0 needs RI-K to be tractable (plain RI-J did not finish 2 SCF
iterations in 4 minutes in testing) -- this script always enables it,
unlike prepare_functional_benchmark.py's optional --enable-rik (whose
hardcoded `define` stdin sequence assumes a fixed number of "delete
leftover data group?" prompts that didn't match this molecule's control
files; this script uses pattern-matching pexpect navigation instead of
counting blank lines, which is robust to that).

Single-point on an already-optimized geometry (no $optimize/$statpt
block), so there is no force/gradient-convergence parameter here -- see
setup_metal_complex_template.py or setup_r2scan3c_template.py for that.

Usage:
    python ~/work_turbomole/scripts/setup_pbe0_svp_template.py \
        --xyz r2scan3c_optimized.xyz --charge -1 --nstates 10 \
        --output-dir template_dir

    # sweeping SCF convergence/grid for benchmarking:
    python ~/work_turbomole/scripts/setup_pbe0_svp_template.py \
        --xyz r2scan3c_optimized.xyz --charge -1 --scfconv 8 --grid m5 \
        --output-dir template_dir_m5
"""

import argparse
import os
import subprocess
import sys
import time


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Build a reusable PBE0/def2-SVP RI-K + TD-DFT "
                    "Turbomole template for one protonation state.")
    parser.add_argument("--xyz", required=True,
                         help="Representative .xyz for this state "
                              "(geometry itself does not matter for the "
                              "template -- only composition/charge do)")
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--nstates", type=int, default=10,
                         help="Number of TD-DFT singlet roots [default: 10]")
    parser.add_argument("--scfconv", type=int, default=7,
                         help="SCF energy convergence criterion, as the N "
                              "in $scfconv N (10^-N Hartree) [default: 7, "
                              "matches define's own default -- unlike the "
                              "HSE06/metal-complex production templates, "
                              "this script did not previously override it].")
    parser.add_argument("--grid", default="m4",
                         help="DFT numerical integration grid size, as "
                              "passed to `grid <value>` in define's dft "
                              "menu (e.g. m3/m4/m5) [default: m4].")
    parser.add_argument("--output-dir", required=True)
    return parser


TDDFT_BLOCK = (
    "$scfconv {scfconv}\n"
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


def run_define(work_dir, charge, grid, log_path):
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
        send("b all def2-SVP", wait=2.0)
        send("*", wait=1.0)
        send("eht")
        send("")
        send(str(charge))
        send("")
        send("dft")
        send("func pbe0")
        send(f"grid {grid}")
        send("on")
        send("")
        send("ri")
        send("on", wait=1.0)
        send("m")
        send("")
        send("")
        advance_to(child, ["GENERAL MENU"])
        send("rijk", wait=1.0)
        advance_to(child, ["ENTER RIJK-OPTION", "jkbas"])
        send("jkbas", wait=2.0)
        send("", wait=1.0)
        send("", wait=1.0)
        send("", wait=1.0)
        send("*", wait=2.0)
        advance_to(child, ["GENERAL MENU"], max_steps=5, step_wait=3)
        send("*", wait=2.0)
        try:
            child.expect(pexpect.EOF, timeout=25)
        finally:
            child.close()


def main():
    args = build_arg_parser().parse_args()

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
    run_define(args.output_dir, args.charge, args.grid, log_path)

    control_path = os.path.join(args.output_dir, "control")
    if not os.path.isfile(control_path):
        raise RuntimeError(f"define did not produce a control file -- check {log_path}")

    with open(control_path) as f:
        text = f.read()
    if "pbe0" not in text or "$rik" not in text:
        raise RuntimeError(
            f"control does not show the expected pbe0/RI-K setup -- "
            f"check {log_path} and {control_path} by hand.")

    text = text.replace(
        "$end",
        TDDFT_BLOCK.format(nstates=args.nstates, scfconv=args.scfconv) + "$end",
        1,
    )
    with open(control_path, "w") as f:
        f.write(text)

    print(f"Template ready in {args.output_dir}/ (control, basis, auxbasis, mos, "
          f"coord; charge={args.charge}, pbe0/def2-SVP, RI-J+RI-K, grid "
          f"{args.grid}, scfconv {args.scfconv}, TD-DFT {args.nstates} "
          f"singlet states).")


if __name__ == "__main__":
    main()
