#!/usr/bin/env python3
"""
setup_metal_complex_template.py

Drives Turbomole's `define` (via pexpect) to build one production-level
Turbomole job directory (control + basis + auxbasis + mos/alpha+beta)
for a metal-dye complex (Cu/Ni/Zn-alizarin-type systems), for direct
geometry/property comparison against an analogous VASP calculation.

Handles both open-shell (UHF, e.g. the Cu(II) doublet -- an odd total
electron count) and closed-shell (RHF, e.g. Ni/Zn singlets) systems
automatically: `define` itself detects the electron-count parity from
the coordinates/charge and switches to UHF with the lowest-multiplicity
aufbau occupation on its own (confirmed empirically for the neutral
Cu complex here: 335 electrons -> "FOUND HALF-OPEN SHELL CASE WITH
MULTIPLICITY 2" without any manual multiplicity input). This script
therefore never types a multiplicity -- it just accepts define's
defaults at every occupation-related prompt via the pattern-matching
`advance_to` helper (GENERAL MENU banner), which is robust regardless
of how many extra "DEFAULT=y" confirmations a transition-metal EHT
guess triggers (Cu needed 2 extra vs. a plain main-group atom; Ni/Zn
were not separately counted -- this is why a fixed blank-count script
would be fragile here, unlike the simpler organic-molecule templates
this is adapted from).

Functional-specific RI treatment (a numerical-acceleration choice, not
a physics/comparability difference -- basis, grid, dispersion, charge,
and SCF convergence stay identical across functionals so the functional
is the only thing that varies, per the user's explicit requirement):
  - hse06: RI-J only (screened exchange via $senex, added by text edit,
    same as Paper 1's artepilin C production template).
  - pbe:   RI-J only (pure GGA, no exact exchange).
  - pbe0:  RI-J + RI-K (global hybrid; plain RI-J does not converge in
           reasonable time -- confirmed previously for PBE0/def2-SVP,
           same mechanism applies here). Uses def2-TZVP (not def2-SVP)
           to match hse06/pbe's basis for comparability -- def2-TZVP is
           a mainstream, fully RI-JK-covered basis, unlike def2-mTZVP.

Usage:
    python ~/work_turbomole/scripts/setup_metal_complex_template.py \
        --coord coord --charge 0 --functional hse06 \
        --nstates 30 --output-dir hse06_m4_d3

    # sweeping SCF/grid/force-convergence for benchmarking:
    python ~/work_turbomole/scripts/setup_metal_complex_template.py \
        --coord coord --charge 0 --functional pbe0 \
        --scfconv 7 --grid m5 --thrmaxgrad 1.0d-3 \
        --output-dir pbe0_m5_d3_loose
"""

import argparse
import os
import re
import shutil
import sys
import time


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Build a production Turbomole template (control/basis/"
                    "auxbasis/mos) for one metal-complex system/functional.")
    parser.add_argument("--coord", required=True,
                         help="Path to an existing Turbomole 'coord' file "
                              "for this system.")
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--functional", required=True,
                         choices=["hse06", "pbe", "pbe0"])
    parser.add_argument("--nstates", type=int, default=30,
                         help="Number of TD-DFT singlet roots (1..N) "
                              "[default: 30]")
    parser.add_argument("--scfconv", type=int, default=8,
                         help="SCF energy convergence criterion, as the N "
                              "in $scfconv N (10^-N Hartree) [default: 8, "
                              "matches Paper 1 production].")
    parser.add_argument("--grid", default="m4",
                         help="DFT numerical integration grid size, as "
                              "passed to `grid <value>` in define's dft "
                              "menu (e.g. m3/m4/m5) [default: m4].")
    parser.add_argument("--thrmaxgrad", default="5.0d-4",
                         help="Geometry optimization force-convergence "
                              "criterion: $statpt's thrmaxgrad (max "
                              "gradient component, Hartree/bohr) "
                              "[default: 5.0d-4, Paper 1 production value]. "
                              "The other four $statpt thresholds "
                              "(threchange, thrrmsgrad, thrmaxdispl, "
                              "thrrmsdispl) stay fixed at Paper 1's values.")
    parser.add_argument("--output-dir", required=True)
    return parser


STATPT_BLOCK = (
    "$statpt\n"
    "   threchange   1.0d-7\n"
    "   thrmaxgrad   {thrmaxgrad}\n"
    "   thrrmsgrad   2.0d-4\n"
    "   thrmaxdispl  5.0d-4\n"
    "   thrrmsdispl  2.0d-4\n"
)

COMMON_EXTRA_BLOCK = (
    "$disp3 bj\n"
    "$scfconv {scfconv}\n"
    "$denconv 1d-7\n"
    "$scfiterlimit 150\n"
    "$scfinstab rpas\n"
    "$soes\n"
    " a           1-{nstates}\n"
    "$rpacor 500\n"
) + STATPT_BLOCK

HSE06_EXTRA = "$senex\n$rij\n" + COMMON_EXTRA_BLOCK
PBE_EXTRA = "$rij\n" + COMMON_EXTRA_BLOCK
PBE0_EXTRA = COMMON_EXTRA_BLOCK  # $rij + $rik come from the `ri`/`rijk` menu steps


def advance_to(child, patterns, max_steps=30, step_wait=3):
    import pexpect
    for _ in range(max_steps):
        idx = child.expect(patterns + [pexpect.TIMEOUT], timeout=step_wait)
        if idx < len(patterns):
            return idx
        child.sendline("")
    raise RuntimeError(f"define: did not reach any of {patterns}")


def run_define(work_dir, charge, functional, grid, log_path):
    import pexpect
    child = pexpect.spawn("define", cwd=work_dir, timeout=25, encoding="utf-8")
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
        send("eht", wait=1.0)
        send(str(charge), wait=1.0)
        # Robust to any number of extra transition-metal EHT-guess /
        # natural-orbital confirmations between charge entry and the
        # occupation summary landing back in GENERAL MENU.
        advance_to(child, ["GENERAL MENU"], max_steps=30, step_wait=4)

        send("dft", wait=1.0)
        send(f"func {functional}", wait=1.0)
        send(f"grid {grid}", wait=1.0)
        send("on", wait=1.5)
        send("", wait=1.0)

        send("ri", wait=1.0)
        send("on", wait=1.5)
        send("m", wait=1.0)
        send("", wait=1.0)
        send("", wait=1.0)

        if functional == "pbe0":
            advance_to(child, ["GENERAL MENU"], max_steps=10, step_wait=3)
            send("rijk", wait=1.0)
            advance_to(child, ["ENTER RIJK-OPTION", "jkbas"], max_steps=10, step_wait=3)
            send("jkbas", wait=2.0)
            send("", wait=1.0)
            send("", wait=1.0)
            send("", wait=1.0)

        send("*", wait=2.0)
        advance_to(child, ["GENERAL MENU"], max_steps=10, step_wait=3)
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

    coord_dest = os.path.join(args.output_dir, "coord")
    shutil.copyfile(args.coord, coord_dest)

    log_path = os.path.join(args.output_dir, "define_session.log")
    run_define(args.output_dir, args.charge, args.functional, args.grid, log_path)

    control_path = os.path.join(args.output_dir, "control")
    if not os.path.isfile(control_path):
        raise RuntimeError(f"define did not produce a control file -- check {log_path}")

    with open(control_path) as f:
        text = f.read()
    if args.functional not in text:
        raise RuntimeError(
            f"control does not show the expected {args.functional} setup -- "
            f"check {log_path} and {control_path} by hand.")
    if args.functional == "pbe0" and "$rik" not in text:
        raise RuntimeError(
            f"control is missing $rik for pbe0 -- check {log_path} and "
            f"{control_path} by hand.")

    # define's own default $scfconv would otherwise duplicate the
    # $scfconv override added below.
    text = re.sub(r"^\$scfconv\s+\d+\s*\n", "", text, count=1, flags=re.MULTILINE)
    # Same for its default $scfiterlimit (30).
    text = re.sub(r"^\$scfiterlimit\s+\d+\s*\n", "", text, count=1, flags=re.MULTILINE)

    extra = {"hse06": HSE06_EXTRA, "pbe": PBE_EXTRA, "pbe0": PBE0_EXTRA}[args.functional]
    text = text.replace(
        "$end",
        extra.format(nstates=args.nstates, scfconv=args.scfconv,
                     thrmaxgrad=args.thrmaxgrad) + "$end",
        1,
    )
    with open(control_path, "w") as f:
        f.write(text)

    print(f"Template ready in {args.output_dir}/ (control, basis, auxbasis, "
          f"mos/alpha+beta, coord; charge={args.charge}, "
          f"{args.functional}-D3(BJ)/def2-TZVP, grid {args.grid}, scfconv "
          f"{args.scfconv}, thrmaxgrad {args.thrmaxgrad}, TD-DFT states "
          f"1-{args.nstates}).")


if __name__ == "__main__":
    main()
