#!/usr/bin/env python3
"""
prepare_functional_benchmark.py

Builds a single-point ridft+escf folder for a TD-DFT functional
benchmark, starting from an already-converged geometry-optimization
folder (coord, basis, auxbasis, mos). Keeps only the data groups
needed for a fixed-geometry SCF + TD-DFT run (RI-J, grid, dispersion,
occupation, TD-DFT roots) and drops everything specific to the
optimizer/Hessian run that produced the source folder (redundant
internal coordinates, $forceupdate, hessian/vibrational data, stale
$last-step / $xctype / $dipole metadata from the source functional).

The source `mos` is copied over as an SCF initial guess (a converged
wavefunction from a nearby functional restarts faster than an EHT
guess); ridft/escf reconverge it for the new functional from scratch.

Global hybrids without screened exchange (PBE0, B3LYP) need RI-K
(--enable-rik) to be tractable: with only RI-J, exact exchange falls
back to conventional 4-index integrals every SCF iteration, which is
dramatically slower than screened functionals like HSE06 (confirmed
empirically -- a plain-RI-J PBE0 single point did not finish one SCF
cycle in 3 minutes on this molecule; HSE06 -proper finished in 35s).
RI-K needs a matching auxiliary basis in $jkbas; the same def2-TZVP
RI-JK set from Turbomole's jkbasen library is used.

Range-separated hybrids (CAM-B3LYP, wB97X, ...) CANNOT use RI-K --
Turbomole refuses outright ("range-separated hybrids not possible
with RI-K!"), confirmed empirically when --enable-rik was tried on
CAM-B3LYP. Use --enable-senex instead: seminumerical exchange, on for
both SCF iterations and the CPHF/CPKS-type response equations that
TD-DFT (escf) solves -- this accelerates both the ground-state SCF
and the excited-state step, unlike RI-K which only speeds up SCF.

Usage:
    python ~/scripts/prepare_functional_benchmark.py \
        --source-dir azo_trans_hse06_m4_d3 \
        --output-dir azo_trans_pbe0_m4_d3 \
        --dft-keyword pbe0 --nstates 10 --enable-rik

    python ~/scripts/prepare_functional_benchmark.py \
        --source-dir azo_trans_hse06_m4_d3 \
        --output-dir azo_trans_camb3lyp_m4_d3 \
        --dft-keyword cam-b3lyp --nstates 10 --enable-senex
"""

import argparse
import os
import re
import shutil
import subprocess

KEEP_PREFIXES = [
    "title", "symmetry", "coord", "atoms", "basis", "scfmo",
    "closed shells", "scfiterlimit", "thize", "thime", "scfdamp",
    "scfintunit", "scfdiis", "maxcor", "scforbitalshift", "energy",
    "dft", "scfconv", "ricore", "rij", "jbas", "disp3",
    "rundimensions", "scfinstab", "soes", "end",
]

# Validated interactively against Turbomole V7-8-1's `define` (see
# module docstring): accept default title -> accept geometry/basis
# defaults -> keep existing MOs -> don't delete leftover data groups
# -> "rijk" menu -> "jkbas" assigns the matching RI-JK aux basis to
# every atom automatically -> "bl" lists it for a sanity check -> "*"
# three times saves and exits. `define` also re-adds its own default
# optimizer/derivative data groups when it doesn't find any in the
# input control -- those are stripped back out by a second filter
# pass after this runs (see enable_rik()).
DEFINE_RIK_STDIN = "\n\n\nn\nn\nrijk\njkbas\nbl\n*\n*\n*\n"
RIK_KEEP_PREFIXES = ["rik", "jkbas"]

# Validated the same way: general menu -> "senex" -> "on" (switch senex
# on for SCF) -> "n" to the gradient follow-up (no geometry/egrad work
# here) -> "y" to also use it in CP-HF/KS (the response equations
# escf's TD-DFT solves) -> blank to leave the senex menu -> "*" to
# save and exit define.
DEFINE_SENEX_STDIN = "\n\n\nn\nn\nsenex\non\nn\ny\n\n*\n"
SENEX_KEEP_PREFIXES = ["senex", "esenex"]


def run_define_and_refilter(output_dir, stdin, extra_keep_prefixes, label):
    """Runs `define` with a validated stdin script to toggle one
    exchange-acceleration option on, then strips back out the
    optimizer/derivative data groups define re-adds by default when it
    doesn't find any in the input control (even from a killed/aborted
    define session earlier in the same folder -- define writes control
    incrementally, not just at the end, so this filter runs
    unconditionally rather than trusting a clean exit)."""
    result = subprocess.run(
        ["define"], cwd=output_dir, input=stdin,
        capture_output=True, text=True, timeout=60,
    )
    combined_output = result.stdout + result.stderr
    with open(os.path.join(output_dir, f"define_{label}.out"), "w") as f:
        f.write(combined_output)

    if "define ended normally" not in combined_output:
        raise RuntimeError(
            f"define did not end normally while enabling {label} in {output_dir} -- "
            f"check {output_dir}/define_{label}.out before trusting this folder.")

    with open(os.path.join(output_dir, "control"), "r") as f:
        control_text = f.read()

    blocks = split_data_groups(control_text)
    all_prefixes = KEEP_PREFIXES + extra_keep_prefixes
    kept, dropped = [], []
    for block in blocks:
        header = block_header(block)
        matched = any(re.match(re.escape(p) + r"(\s|$)", header) for p in all_prefixes)
        (kept if matched else dropped).append(block)

    print(f"define enabled {label}; re-filtered to {len(kept)} data groups, "
          f"dropped {len(dropped)} that define re-added on its own:")
    for block in dropped:
        print(f"  dropped: ${block_header(block).splitlines()[0]}")

    with open(os.path.join(output_dir, "control"), "w") as f:
        f.write("".join(kept))


def split_data_groups(control_text):
    """Splits a control file into (header, block_text) pairs, one per
    top-level '$...' data group (continuation lines included)."""
    blocks = []
    current_lines = []
    for line in control_text.splitlines(keepends=True):
        if line.startswith("$"):
            if current_lines:
                blocks.append("".join(current_lines))
            current_lines = [line]
        else:
            current_lines.append(line)
    if current_lines:
        blocks.append("".join(current_lines))
    return blocks


def block_header(block):
    first_line = block.splitlines()[0]
    return first_line[1:].strip()  # strip leading '$'


def filter_and_edit_blocks(blocks, dft_keyword, nstates):
    kept, dropped = [], []
    for block in blocks:
        header = block_header(block)
        matched = any(re.match(re.escape(p) + r"(\s|$)", header) for p in KEEP_PREFIXES)
        if not matched:
            dropped.append(header)
            continue

        if header.startswith("dft"):
            block = re.sub(r"(functional\s+)\S+", rf"\g<1>{dft_keyword}", block)
        elif header.startswith("soes"):
            block = re.sub(r"(\ba\s+)\d+", rf"\g<1>{nstates}", block)

        kept.append(block)
    return kept, dropped


def main():
    parser = argparse.ArgumentParser(
        description="Prepare a single-point ridft+escf folder for a TD-DFT "
                    "functional benchmark, from an already-optimized geometry.")
    parser.add_argument("--source-dir", required=True,
                         help="Folder with the converged geometry (control, coord, "
                              "basis, auxbasis, mos)")
    parser.add_argument("--output-dir", required=True,
                         help="New folder to create for this functional's run")
    parser.add_argument("--dft-keyword", required=True,
                         help="Turbomole $dft functional keyword, e.g. pbe0, b3-lyp, cam-b3lyp")
    parser.add_argument("--nstates", type=int, default=10,
                         help="Number of singlet roots for $soes")
    parser.add_argument("--enable-rik", action="store_true",
                         help="Run `define` to add $rik + a matching RI-JK auxiliary "
                              "basis (needed for global hybrids without screened "
                              "exchange, e.g. PBE0, B3LYP -- NOT usable with "
                              "range-separated hybrids like CAM-B3LYP, Turbomole "
                              "refuses those; use --enable-senex for those instead)")
    parser.add_argument("--enable-senex", action="store_true",
                         help="Run `define` to switch on seminumerical exchange for "
                              "SCF and CPHF/CPKS (needed for range-separated hybrids "
                              "like CAM-B3LYP, which reject --enable-rik)")
    args = parser.parse_args()

    if args.enable_rik and args.enable_senex:
        raise ValueError("--enable-rik and --enable-senex are alternatives for "
                          "different functional families -- pick one.")

    if os.path.exists(args.output_dir):
        raise FileExistsError(f"{args.output_dir} already exists -- refusing to overwrite.")

    with open(os.path.join(args.source_dir, "control"), "r") as f:
        control_text = f.read()

    blocks = split_data_groups(control_text)
    kept, dropped = filter_and_edit_blocks(blocks, args.dft_keyword, args.nstates)

    print(f"Kept {len(kept)} data groups, dropped {len(dropped)}:")
    for header in dropped:
        print(f"  dropped: ${header.splitlines()[0] if header else header}")

    os.makedirs(args.output_dir)
    with open(os.path.join(args.output_dir, "control"), "w") as f:
        f.write("".join(kept))

    for fname in ["coord", "basis", "auxbasis", "mos"]:
        shutil.copy(os.path.join(args.source_dir, fname),
                    os.path.join(args.output_dir, fname))

    if args.enable_rik:
        run_define_and_refilter(args.output_dir, DEFINE_RIK_STDIN, RIK_KEEP_PREFIXES, "rik")
    elif args.enable_senex:
        run_define_and_refilter(args.output_dir, DEFINE_SENEX_STDIN, SENEX_KEEP_PREFIXES, "senex")

    print(f"\nNew folder ready -> {args.output_dir}")
    print("(control + coord + basis + auxbasis + mos; mos is a warm-start guess "
          "from the source functional, ridft will reconverge it)")


if __name__ == "__main__":
    main()
