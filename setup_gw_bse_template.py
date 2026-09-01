#!/usr/bin/env python3
"""
setup_gw_bse_template.py

Drives Turbomole's `define` (via pexpect) to build a two-stage RI-GW +
GW-BSE/TDA template: one G0W0 quasiparticle-energy calculation
(directory `gw/`) and one Bethe-Salpeter exciton calculation on top of
it (directory `bse/`), for a single molecule/charge/functional.

This fills in a gap flagged in MANUAL.md's Estagio 1d ("o procedimento
exato dessa rodada (modulo GW do Turbomole) ainda nao esta documentado
neste repositorio") -- the downstream parsers/plots
(parse_qpenergies_dos.py, compute_exciton_binding.py, compare_low_states.py)
already existed and already assumed `qpenergies.dat` + `escf_bse.out`,
but nothing in this repo actually produced them until this script.

Discovered empirically (via `define`'s GENERAL MENU -> `gw` submenu,
confirmed by running a real (tiny, water) test end-to-end -- no
existing example of this workflow existed anywhere in this repo to
copy from):

  1. GW and BSE are two SEPARATE `escf` runs with two DIFFERENT
     control files, not one combined run:
       - `gw/control`: ground-state dft/ri setup + `$rigw` (RI-accelerated
         G0W0, analytic continuation, N^4 -- NOT the plain `$gw` spectral-
         function option, which is N^6 and was not tested here) with
         `fullspec` (quasiparticle energies for ALL orbitals, not just
         near the gap) -- `ridft` then `escf` here produces
         `qpenergies.dat`.
       - `bse/control`: the SAME ground-state dft/ri setup (no `$rigw`
         block at all) plus `$scfinstab rpas` + `$bse` + `$soes` (number
         of excitonic states). This step needs `qpenergies.dat` copied
         into its directory (from the finished `gw/` run) before
         `ridft`+`escf` -- confirmed by inspecting escf_bse.out's
         "Dominant contributions" orbital energies: they exactly matched
         qpenergies.dat's G0W0-corrected values (e.g. -18.28 eV / 4.52 eV
         for water's HOMO-1/LUMO), NOT the bare Kohn-Sham ones
         (-12.53 eV / 0.81 eV) -- i.e. BSE really is reading and using
         the GW correction, not silently falling back to plain TD-DFT.
  2. define's `gw` submenu quirk: bare `rigw` does NOT toggle it on
     (silently re-displays the same "off" menu, no error) -- it must be
     `rigw on` explicitly. This is inconsistent with the plain `dft`/`ri`
     submenus (where `on` is a separate follow-up line) and with the
     `ex` menu's `bse` toggle (which DOES work as a bare command) --
     confirmed by testing all three side by side; don't assume one
     on/off convention generalizes to every define submenu.
  3. `fullspec` removes GW's default `ips+`/`gap` restriction (per its
     own menu text: "CALCULATE ALL QP ENERGIES (REMOVES ips+ AND gap)"),
     giving a full orbital-by-orbital qpenergies.dat table -- required
     for parse_qpenergies_dos.py's DOS comparison use case (MANUAL.md's
     stated goal), not just a HOMO-LUMO gap number. Costlier for large
     systems (scales with total orbital count); use --no-fullspec for
     just the frontier-orbital region on bigger molecules.
  4. RI-GW needs RI-J (already on for any of this project's functionals)
     and pulls in `$rick` automatically -- no extra RI-JK setup command
     was needed for the `gw` step itself, unlike PBE0's OWN ground-state
     RI-K requirement (see setup_metal_complex_template.py) which is a
     separate, independent RI-K use.

NOT yet tested here (do a small test before trusting on a real system):
  - PBE0 as the mean-field starting point together with RI-GW (only
    plain PBE and HSE06 were run end-to-end during this exploration --
    PBE0 additionally needs its own ground-state RI-K/rijk menu step,
    reused from setup_metal_complex_template.py, but the RI-GW + RI-K
    combination itself was not exercised).
  - Open-shell (UHF) systems -- GW/BSE literature and this Turbomole
    module are predominantly closed-shell; this script assumes a
    closed-shell molecule (even electron count) and does not attempt
    setup_metal_complex_template.py's UHF-robust `advance_to` occupation
    handling.
  - Anything beyond a 3-atom test molecule -- cost/scaling on a real
    system (e.g. artepilin C, ~600 AOs) is unknown; benchmark before
    committing cluster walltime.

Usage:
    python ~/work_turbomole/scripts/setup_gw_bse_template.py \
        --xyz structure.xyz --charge 0 --functional hse06 \
        --nstates 10 --output-dir gwbse_test

    # then, on the cluster:
    cd gwbse_test/gw   && ridft > ridft.out && escf > escf.out
    cp qpenergies.dat ../bse/
    cd ../bse && ridft > ridft.out && escf > escf_bse.out
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import time


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Build a two-stage RI-GW + GW-BSE/TDA Turbomole "
                    "template (gw/ + bse/ subdirectories) for one "
                    "molecule/charge/functional.")
    parser.add_argument("--xyz", required=True,
                         help="Geometry for the mean-field (DFT) reference "
                              "-- an already-optimized structure, since "
                              "this builds single-point property "
                              "calculations only (no $optimize/$statpt).")
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--functional", default="hse06",
                         choices=["hse06", "pbe", "pbe0"],
                         help="Mean-field starting point for G0W0 "
                              "[default: hse06]. hse06/pbe confirmed "
                              "working with RI-GW; pbe0 (RI-K ground "
                              "state) has NOT been tested together with "
                              "RI-GW -- verify on a small system first.")
    parser.add_argument("--nstates", type=int, default=10,
                         help="Number of BSE/TDA excitonic singlet states "
                              "(the $soes count in bse/control) "
                              "[default: 10].")
    parser.add_argument("--scfconv", type=int, default=8,
                         help="SCF energy convergence criterion, as the N "
                              "in $scfconv N (10^-N Hartree) [default: 8, "
                              "matches Paper 1 production].")
    parser.add_argument("--grid", default="m4",
                         help="DFT numerical integration grid size, as "
                              "passed to `grid <value>` in define's dft "
                              "menu (e.g. m3/m4/m5) [default: m4].")
    parser.add_argument("--no-fullspec", dest="fullspec", action="store_false",
                         help="Restrict GW to the default near-gap orbital "
                              "window (ips+) instead of ALL orbitals -- "
                              "cheaper on large systems, but insufficient "
                              "for a full DOS comparison (MANUAL.md "
                              "Estagio 1d's use case).")
    parser.add_argument("--output-dir", required=True)
    return parser


GW_EXTRA_BLOCK = (
    "$disp3 bj\n"
    "$scfconv {scfconv}\n"
    "$denconv 1d-7\n"
)

BSE_EXTRA_BLOCK = (
    "$disp3 bj\n"
    "$scfconv {scfconv}\n"
    "$denconv 1d-7\n"
    "$scfinstab rpas\n"
    "$bse\n"
    "$soes\n"
    " a           {nstates}\n"
    "$rpacor 500\n"
)

HSE06_GROUND_STATE = "$senex\n$rij\n"
PBE_GROUND_STATE = "$rij\n"
PBE0_GROUND_STATE = ""  # $rij + $rik come from the `ri`/`rijk` menu steps


def advance_to(child, patterns, max_steps=30, step_wait=3):
    import pexpect
    for _ in range(max_steps):
        idx = child.expect(patterns + [pexpect.TIMEOUT], timeout=step_wait)
        if idx < len(patterns):
            return idx
        child.sendline("")
    raise RuntimeError(f"define: did not reach any of {patterns}")


def run_ground_state_steps(send, charge, functional, grid):
    """Geometry/basis/occupation/dft/ri steps shared by both gw/ and bse/."""
    send("")
    send("")
    send("a coord")
    send("*")
    send("no")
    send("b all def2-TZVP", wait=2.0)
    send("*", wait=1.0)
    send("eht", wait=1.0)
    send(str(charge), wait=1.0)
    send("", wait=1.5)
    send("", wait=1.5)

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


def build_gw_dir(work_dir, charge, functional, grid, fullspec, scfconv, log_path):
    import pexpect
    child = pexpect.spawn("define", cwd=work_dir, timeout=25, encoding="utf-8")
    with open(log_path, "w") as logf:
        child.logfile = logf

        def send(s, wait=0.8):
            child.sendline(s)
            time.sleep(wait)

        run_ground_state_steps(send, charge, functional, grid)

        if functional == "pbe0":
            advance_to(child, ["GENERAL MENU"], max_steps=10, step_wait=3)
            send("rijk", wait=1.0)
            advance_to(child, ["ENTER RIJK-OPTION", "jkbas"], max_steps=10, step_wait=3)
            send("jkbas", wait=2.0)
            send("", wait=1.0)
            send("", wait=1.0)
            send("", wait=1.0)

        advance_to(child, ["GENERAL MENU"], max_steps=10, step_wait=3)
        send("gw", wait=3.0)
        send("rigw on", wait=5.0)
        if fullspec:
            send("fullspec", wait=2.0)
        send("", wait=2.0)
        send("*", wait=2.0)

        advance_to(child, ["GENERAL MENU"], max_steps=10, step_wait=3)
        send("*", wait=2.0)
        try:
            child.expect(pexpect.EOF, timeout=25)
        finally:
            child.close()

    control_path = os.path.join(work_dir, "control")
    with open(control_path) as f:
        text = f.read()
    if "$rigw" not in text:
        raise RuntimeError(
            f"gw/control is missing $rigw -- check {log_path} and "
            f"{control_path} by hand.")

    ground_state_block = {
        "hse06": HSE06_GROUND_STATE,
        "pbe": PBE_GROUND_STATE,
        "pbe0": PBE0_GROUND_STATE,
    }[functional]
    text = re.sub(r"^\$scfconv\s+\d+\s*\n", "", text, count=1, flags=re.MULTILINE)
    text = text.replace(
        "$end",
        ground_state_block + GW_EXTRA_BLOCK.format(scfconv=scfconv) + "$end",
        1,
    )
    with open(control_path, "w") as f:
        f.write(text)


def build_bse_dir(work_dir, charge, functional, grid, nstates, scfconv, log_path):
    import pexpect
    child = pexpect.spawn("define", cwd=work_dir, timeout=25, encoding="utf-8")
    with open(log_path, "w") as logf:
        child.logfile = logf

        def send(s, wait=0.8):
            child.sendline(s)
            time.sleep(wait)

        run_ground_state_steps(send, charge, functional, grid)

        if functional == "pbe0":
            advance_to(child, ["GENERAL MENU"], max_steps=10, step_wait=3)
            send("rijk", wait=1.0)
            advance_to(child, ["ENTER RIJK-OPTION", "jkbas"], max_steps=10, step_wait=3)
            send("jkbas", wait=2.0)
            send("", wait=1.0)
            send("", wait=1.0)
            send("", wait=1.0)

        advance_to(child, ["GENERAL MENU"], max_steps=10, step_wait=3)
        send("ex", wait=2.0)
        send("bse", wait=3.0)
        send("*", wait=2.0)

        advance_to(child, ["GENERAL MENU"], max_steps=10, step_wait=3)
        send("*", wait=2.0)
        try:
            child.expect(pexpect.EOF, timeout=25)
        finally:
            child.close()

    control_path = os.path.join(work_dir, "control")
    with open(control_path) as f:
        text = f.read()
    if "$bse" not in text:
        raise RuntimeError(
            f"bse/control is missing $bse -- check {log_path} and "
            f"{control_path} by hand.")

    ground_state_block = {
        "hse06": HSE06_GROUND_STATE,
        "pbe": PBE_GROUND_STATE,
        "pbe0": PBE0_GROUND_STATE,
    }[functional]
    text = re.sub(r"^\$scfconv\s+\d+\s*\n", "", text, count=1, flags=re.MULTILINE)
    text = text.replace(
        "$end",
        ground_state_block + BSE_EXTRA_BLOCK.format(scfconv=scfconv, nstates=nstates) + "$end",
        1,
    )
    with open(control_path, "w") as f:
        f.write(text)


def main():
    args = build_arg_parser().parse_args()

    try:
        import pexpect  # noqa: F401
    except ImportError:
        sys.exit("This script requires pexpect: pip install pexpect")

    if os.path.exists(args.output_dir):
        raise FileExistsError(f"{args.output_dir} already exists -- refusing to overwrite.")
    os.makedirs(args.output_dir)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    gw_dir = os.path.join(args.output_dir, "gw")
    bse_dir = os.path.join(args.output_dir, "bse")
    os.makedirs(gw_dir)
    os.makedirs(bse_dir)

    for d in (gw_dir, bse_dir):
        subprocess.run(
            ["python3", os.path.join(script_dir, "xyz_to_coord.py"),
             "--xyz", os.path.abspath(args.xyz),
             "--output", os.path.join(d, "coord")],
            check=True,
        )

    build_gw_dir(gw_dir, args.charge, args.functional, args.grid,
                 args.fullspec, args.scfconv,
                 os.path.join(gw_dir, "define_session.log"))
    build_bse_dir(bse_dir, args.charge, args.functional, args.grid,
                  args.nstates, args.scfconv,
                  os.path.join(bse_dir, "define_session.log"))

    print(f"Template ready in {args.output_dir}/:\n"
          f"  gw/  -- {args.functional}-D3(BJ)/def2-TZVP + RI-GW"
          f"{' (fullspec)' if args.fullspec else ' (ips+ near-gap only)'}"
          f", grid {args.grid}, scfconv {args.scfconv}.\n"
          f"         Run: cd {gw_dir} && ridft > ridft.out && escf > escf.out\n"
          f"         (produces qpenergies.dat)\n"
          f"  bse/ -- same ground state, GW-BSE/TDA, {args.nstates} singlet "
          f"excitonic states.\n"
          f"         Run AFTER gw/ finishes: cp {gw_dir}/qpenergies.dat {bse_dir}/ "
          f"&& cd {bse_dir} && ridft > ridft.out && escf > escf_bse.out\n"
          f"Charge={args.charge}. See this script's docstring for what has "
          f"and has not been tested (pbe0+RI-GW and open-shell systems: NOT "
          f"tested).")


if __name__ == "__main__":
    main()
