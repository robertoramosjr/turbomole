#!/usr/bin/env python3
"""
batch_stda_s1_character_matched.py

Replaces batch_stda_s1.py's "take sTDA's lowest-energy state" rule with
orbital-CHARACTER matching, following the diagnosis in
protonation_and_conformer_investigation.md (oscillator-strength
sections): sTDA robustly places a spurious dark state below the true
bright S1-like transition for the neutral and monoanion_b protonation
states, across the whole ax=0.0-0.30 range, so neither "lowest energy"
nor "largest oscillator strength" (contrast against the next-brightest
state is only ~1.3-5x, never a clean separation) is a safe selection
rule on its own. Orbital character IS reliable: matched HSE06's own
dominant occupied->virtual transition (from escf.out's "Dominant
contributions" table, using its "number of occupied orbitals" line to
identify its own HOMO index) in 18/20 neutral and 11/11 monoanion_b
conformers, restricted to the first --window-states sTDA states by
energy (an unbounded search lets a spurious high-energy state with
large f and wrong character hijack the pick -- see
analyze_stda_state_selection.py's diagnostic).

Per conformer: runs tm2molden + std2 once (same pipeline as
batch_stda_s1.py), computes the local HOMO index from std2's own
"oMOs in TDA" line, and picks the state (within the energy window)
whose contributions include the SAME (occ - HOMO, virt - HOMO) offset
as HSE06's own dominant contribution, at |coefficient| >=
--min-match-coeff. Falls back to the largest-f state in the window
(flagged) when no state matches that offset at all -- this happens for
conformers where HSE06's own S1 is not cleanly HOMO->LUMO either, so
there is no single offset to look for.

IMPORTANT (found 2026-09-11, Sec. 9.25 of
protonation_and_conformer_investigation.md): when the target character
is FRAGMENTED across several near-degenerate sTDA states (common in
this molecule -- roughly half of the neutral state's labeled
conformers showed it), picking the state with the single LARGEST
matching coefficient is not reliable -- two or more states can have
near-tied coefficients (e.g. 0.54 vs 0.53) for utterly different
oscillator strengths (e.g. f=0.06 vs f=0.31, a 5x difference), and the
largest-coefficient tie-break has no reason to land on the physically
brighter one. Among states clearing --min-match-coeff, this version
picks the one with the LARGEST OSCILLATOR STRENGTH, not the largest
matching coefficient -- directly optimizing for the quantity B2 needs
(how bright this character actually is), instead of an arbitrary
coefficient tie-break that does not correlate with brightness.

Usage:
    python ~/work_turbomole/scripts/batch_stda_s1_character_matched.py \
        --job-list job_dirs.txt \
        --r2scan3c-jobs layer_a/conformer_jobs \
        --labeled-jobs layer_b/labeled_jobs \
        --ax 0.25 \
        --window-states 10 \
        --output layer_b/s1_stda_character_matched.csv
"""

import argparse
import csv
import os
import re
import subprocess
import tempfile

STATE_ROW_RE = re.compile(
    r"^\s*(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+(.*)$"
)
CONTRIB_RE = re.compile(r"([-\d.]+)\(\s*(\d+)->\s*(\d+)\)")
OMOS_RE = re.compile(r"oMOs in TDA:\s*(\d+)")
HSE06_OCC_COUNT_RE = re.compile(r"number of occupied orbitals\s*:\s*(\d+)")
HSE06_DOMINANT_RE = re.compile(
    r"^\s*(\d+)\s+a\s+[-\d.]+\s+(\d+)\s+a\s+[-\d.]+\s+([\d.]+)\s*$"
)


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
        states.append({
            "state": int(state), "energy_eV": float(ev), "osc": float(fl),
            "contribs": contribs,
        })
    return states or None


def parse_local_homo_index(std2_stdout):
    m = OMOS_RE.search(std2_stdout)
    return int(m.group(1)) if m else None


def parse_hse06_s1(escf_path):
    homo_index = None
    with open(escf_path) as f:
        text = f.read()
    m = HSE06_OCC_COUNT_RE.search(text)
    if m:
        homo_index = int(m.group(1))

    lines = text.splitlines()
    dom_occ = dom_virt = None
    for i, line in enumerate(lines):
        if "1 singlet a excitation" in line:
            for j in range(i, min(i + 60, len(lines))):
                if "Dominant contributions" in lines[j]:
                    for k in range(j, min(j + 6, len(lines))):
                        dm = HSE06_DOMINANT_RE.match(lines[k])
                        if dm:
                            dom_occ, dom_virt = int(dm.group(1)), int(dm.group(2))
                            break
                    break
            break
    return homo_index, dom_occ, dom_virt


def main():
    parser = argparse.ArgumentParser(
        description="Extract sTDA's S1 proxy by matching HSE06's own "
                    "dominant occupied->virtual orbital offset, instead of "
                    "assuming lowest energy or largest oscillator strength.")
    parser.add_argument("--job-list", required=True)
    parser.add_argument("--r2scan3c-jobs", required=True)
    parser.add_argument("--labeled-jobs", required=True,
                         help="Folder of job_NNNN/escf.out (HSE06 production labels)")
    parser.add_argument("--energy-threshold", type=float, default=8.0)
    parser.add_argument("--ax", type=float, default=0.25)
    parser.add_argument("--window-states", type=int, default=10,
                         help="Restrict character/intensity matching to the "
                              "first N sTDA states by energy [default: 10]")
    parser.add_argument("--min-match-coeff", type=float, default=0.3,
                         help="Minimum |coefficient| for a state's "
                              "contribution to count as matching HSE06's "
                              "target offset [default: 0.3 -- distinguishes "
                              "a genuine, non-negligible admixture from "
                              "small numerical mixing noise]")
    parser.add_argument("--std2-bin", default=os.path.expanduser("~/local_bin/std2"))
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.job_list) as f:
        job_dirs = [line.strip() for line in f if line.strip()]

    rows = []
    for jd in job_dirs:
        src = os.path.join(args.r2scan3c_jobs, jd)
        escf_path = os.path.join(args.labeled_jobs, jd, "escf.out")
        if not os.path.isdir(src) or not os.path.isfile(escf_path):
            print(f"WARNING: missing input for {jd}, skipping")
            continue

        hse06_homo, hse06_dom_occ, hse06_dom_virt = parse_hse06_s1(escf_path)
        if hse06_homo is None or hse06_dom_occ is None:
            print(f"WARNING: could not parse HSE06 reference for {jd}, skipping")
            continue
        occ_offset = hse06_dom_occ - hse06_homo
        virt_offset = hse06_dom_virt - (hse06_homo + 1)

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

        matches = []
        for s in window:
            for coeff, occ_l, virt_l in s["contribs"]:
                if (occ_l - local_homo == occ_offset
                        and virt_l - (local_homo + 1) == virt_offset
                        and abs(coeff) >= args.min_match_coeff):
                    matches.append(s)
                    break

        if matches:
            chosen = max(matches, key=lambda s: s["osc"])
            method = ("character_match" if len(matches) == 1
                       else "character_match_multi_candidate")
        else:
            chosen = max(window, key=lambda s: s["osc"])
            method = "fallback_max_f"

        rows.append({
            "job_dir": jd, "s1_energy_eV_stda": chosen["energy_eV"],
            "s1_osc_stda": chosen["osc"], "selection_method": method,
        })
        print(f"{jd}: E={chosen['energy_eV']:.4f} eV, f={chosen['osc']:.5f} "
              f"[{method}]")

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["job_dir", "s1_energy_eV_stda", "s1_osc_stda",
                           "selection_method"])
        writer.writeheader()
        writer.writerows(rows)

    n_fallback = sum(1 for r in rows if r["selection_method"] == "fallback_max_f")
    print(f"\nExtracted {len(rows)} / {len(job_dirs)} conformers "
          f"({n_fallback} used the max-f fallback).")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
