#!/usr/bin/env python3
"""
check_uncertain_conformer_character.py

For conformers flagged 'uncertain_selection' by
characterize_stda_selection_confidence.py (low max-f/second-f contrast
in the sTDA window, no HSE06 reference available to check orbital
character directly), tests whether the SAME failure signature already
confirmed in the labeled training data (analyze_stda_state_selection.py,
job_0011/job_0072 of the neutral state: the max-f state's dominant
transition is (HOMO-2)->LUMO -- offset (-2, 0) relative to the local
HOMO/LUMO pair -- while a weaker but still competitive state nearby is
the TRUE (HOMO->LUMO, offset (0, 0)) transition -- is present.

This cannot prove a conformer's selection is wrong without a real
HSE06 label, but it can separate two very different situations among
the flagged conformers:
  - the chosen (max-f) state's own dominant transition is ALREADY
    (0, 0) HOMO->LUMO -- the low contrast is against some OTHER,
    chemically different state, and does not indicate the same
    failure mode seen in training; likely a correct pick with just a
    thin margin.
  - the chosen state's dominant transition matches the KNOWN mismatch
    offset (typically (-2, 0), but reported for whatever offset the
    training mismatches actually showed) AND a (0, 0) HOMO->LUMO state
    exists nearby in the window -- the same failure signature as the
    2 known training mismatches; flag as a probable real error.
  - neither pattern -- reported for human judgement.

Usage:
    python ~/work_turbomole/scripts/check_uncertain_conformer_character.py \
        --job-list uncertain_job_dirs.txt \
        --r2scan3c-jobs layer_a/conformer_jobs \
        --known-mismatch-offsets -2,0 \
        --ax 0.25 --window-states 10 \
        --output layer_b/uncertain_character_check.dat
"""

import argparse
import os
import re
import subprocess
import tempfile

STATE_ROW_RE = re.compile(
    r"^\s*(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+(.*)$"
)
CONTRIB_RE = re.compile(r"([-\d.]+)\(\s*(\d+)->\s*(\d+)\)")
OMOS_RE = re.compile(r"oMOs in TDA:\s*(\d+)")


def run_tm2molden(work_dir):
    result = subprocess.run(
        ["tm2molden"], cwd=work_dir, input="\n\n", capture_output=True,
        text=True, timeout=30,
    )
    return "ended normally" in result.stdout or "ended normally" in result.stderr


def run_std2(work_dir, energy_threshold, ax, std2_bin):
    result = subprocess.run(
        [std2_bin, "-f", "molden.input", "-ax", str(ax), "-e", str(energy_threshold)],
        cwd=work_dir, capture_output=True, text=True, timeout=120,
    )
    return result.stdout


def parse_all_states(std2_stdout):
    lines = std2_stdout.splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("state") and "eV" in line and "fL" in line:
            header_idx = i
            break
    if header_idx is None:
        return None
    states = []
    for line in lines[header_idx + 1:]:
        m = STATE_ROW_RE.match(line)
        if not m:
            break
        state, ev, nm, fl, rv, contribs_str = m.groups()
        contribs = [(float(c), int(o), int(v))
                    for c, o, v in CONTRIB_RE.findall(contribs_str)]
        states.append({"state": int(state), "energy_eV": float(ev),
                        "osc": float(fl), "contribs": contribs})
    return states or None


def parse_local_homo_index(std2_stdout):
    m = OMOS_RE.search(std2_stdout)
    return int(m.group(1)) if m else None


def dominant_offset(state, local_homo):
    if not state["contribs"]:
        return None
    coeff, occ_l, virt_l = max(state["contribs"], key=lambda c: abs(c[0]))
    return (occ_l - local_homo, virt_l - (local_homo + 1)), coeff


def main():
    parser = argparse.ArgumentParser(
        description="Check whether uncertain_selection conformers show the "
                    "same orbital-offset failure signature as the known "
                    "training mismatches, without needing an HSE06 label.")
    parser.add_argument("--job-list", required=True)
    parser.add_argument("--r2scan3c-jobs", required=True)
    parser.add_argument("--known-mismatch-offsets", default="-2,0",
                         help="occ_offset,virt_offset of the known training "
                              "mismatch's WRONGLY-chosen state, comma-"
                              "separated [default: -2,0, from job_0011/"
                              "job_0072 of the neutral state]")
    parser.add_argument("--energy-threshold", type=float, default=8.0)
    parser.add_argument("--ax", type=float, default=0.25)
    parser.add_argument("--window-states", type=int, default=10)
    parser.add_argument("--top-k", type=int, default=3,
                         help="Report the top-K states by oscillator "
                              "strength within the window [default: 3]")
    parser.add_argument("--std2-bin", default=os.path.expanduser("~/local_bin/std2"))
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    known_occ_off, known_virt_off = (int(x) for x in args.known_mismatch_offsets.split(","))

    with open(args.job_list) as f:
        job_dirs = [line.strip() for line in f if line.strip()]

    lines_out = []
    for jd in job_dirs:
        src = os.path.join(args.r2scan3c_jobs, jd)
        if not os.path.isdir(src):
            print(f"WARNING: {src} not found, skipping")
            continue
        with tempfile.TemporaryDirectory() as tmp:
            for fname in ("control", "basis", "coord", "mos"):
                src_file = os.path.join(src, fname)
                if os.path.isfile(src_file):
                    subprocess.run(["cp", src_file, tmp], check=True)
            if not run_tm2molden(tmp):
                print(f"WARNING: tm2molden failed for {jd}")
                continue
            stdout = run_std2(tmp, args.energy_threshold, args.ax, args.std2_bin)

        states = parse_all_states(stdout)
        local_homo = parse_local_homo_index(stdout)
        if states is None or local_homo is None:
            print(f"WARNING: could not parse std2 output for {jd}")
            continue

        window = sorted(states, key=lambda s: s["energy_eV"])[:args.window_states]
        top = sorted(window, key=lambda s: s["osc"], reverse=True)[:args.top_k]

        chosen = top[0]
        chosen_offset, chosen_coeff = dominant_offset(chosen, local_homo)
        is_chosen_homo_lumo = chosen_offset == (0, 0)
        is_chosen_known_mismatch = chosen_offset == (known_occ_off, known_virt_off)

        has_homo_lumo_competitor = False
        for s in top[1:]:
            off, _ = dominant_offset(s, local_homo)
            if off == (0, 0):
                has_homo_lumo_competitor = True

        if is_chosen_homo_lumo:
            verdict = "LIKELY_CORRECT (chosen state is already HOMO->LUMO)"
        elif is_chosen_known_mismatch and has_homo_lumo_competitor:
            verdict = "PROBABLE_MISMATCH (matches known training failure signature)"
        else:
            verdict = "UNCLEAR (neither known pattern -- needs human judgement)"

        detail = "; ".join(
            f"f={s['osc']:.4f} E={s['energy_eV']:.3f}eV offset={dominant_offset(s, local_homo)[0]}"
            for s in top
        )
        line = f"{jd}: chosen_offset={chosen_offset} [{verdict}]  top-{args.top_k}: {detail}"
        print(line)
        lines_out.append(line)

    with open(args.output, "w") as f:
        f.write("\n".join(lines_out) + "\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
