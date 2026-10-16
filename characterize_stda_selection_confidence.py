#!/usr/bin/env python3
"""
characterize_stda_selection_confidence.py

Production-time confidence check for B2's sTDA oscillator-strength
proxy, for conformers that have NO HSE06 reference (the full Layer A
ensemble, beyond the small labeled subset used to close B1/B2's
feature selection). No HSE06 means no orbital-character validation is
possible (that needs HSE06's own dominant occupied->virtual pair) --
the only sTDA-only signal available is the max-f/second-max-f CONTRAST
within the energy window (see analyze_stda_state_selection.py), which
Task 2's diagnostic showed is only mildly discriminating (median
1.9-5.3x, never above ~10x) but still correlates with whether the
max-f pick's orbital character actually matched HSE06 in the labeled
set.

IMPORTANT SCOPE NOTE: this max-f criterion and its contrast-based
confidence check are only meaningful for protonation states where the
true S1 is BRIGHT (neutral, monoanion_b) -- for monoanion_a and
dianion, HSE06's own S1 is already near-dark (f~0.0003-0.0006), so
"largest f in the window" is selecting the WRONG transition by
construction; those two states never had a state-selection problem in
the first place (their plain "state 1 by energy rank" already agreed
well with HSE06 -- see the original oscillator-strength diagnostic).
Do NOT run this on monoanion_a/dianion conformers; it would produce a
confidence score for a criterion that does not apply there.

Given --calibration-csv (analyze_stda_state_selection.py's
--max-window-states output, which has both the contrast ratio and
whether the pick's orbital character matched HSE06 for the labeled
subset), derives a threshold as the highest contrast ratio among the
known CHARACTER_MISMATCH cases -- any production conformer at or below
that contrast is flagged "uncertain": its character has an empirically
demonstrated chance of being wrong, on the labeled evidence.

Usage:
    python ~/work_turbomole/scripts/characterize_stda_selection_confidence.py \
        --job-list job_dirs.txt \
        --r2scan3c-jobs layer_a/conformer_jobs \
        --calibration-csv layer_b/stda_state_selection_check_win10.csv \
        --ax 0.25 --window-states 10 \
        --output layer_b/stda_confidence_full_ensemble.csv
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
        state, ev, nm, fl, rv, _ = m.groups()
        states.append({"energy_eV": float(ev), "osc": float(fl)})
    return states or None


def derive_threshold(calibration_csv):
    mismatched_contrasts = []
    matched_contrasts = []
    with open(calibration_csv) as f:
        for row in csv.DictReader(f):
            ratio_str = row["stda_f_contrast_ratio"]
            if ratio_str == "inf":
                continue
            ratio = float(ratio_str)
            if "CHARACTER_MISMATCH" in row["flags"]:
                mismatched_contrasts.append(ratio)
            else:
                matched_contrasts.append(ratio)
    if not mismatched_contrasts:
        return None, matched_contrasts, mismatched_contrasts
    threshold = max(mismatched_contrasts)
    return threshold, matched_contrasts, mismatched_contrasts


def main():
    parser = argparse.ArgumentParser(
        description="Compute sTDA max-f/window contrast for every conformer "
                    "in the full (unlabeled) Layer A ensemble, and flag "
                    "conformers whose contrast falls at or below the "
                    "highest contrast seen among known character mismatches "
                    "in the labeled calibration set.")
    parser.add_argument("--job-list", required=True)
    parser.add_argument("--r2scan3c-jobs", required=True)
    parser.add_argument("--calibration-csv", required=True,
                         help="analyze_stda_state_selection.py output for "
                              "the labeled subset of the SAME protonation "
                              "state (must include stda_f_contrast_ratio "
                              "and flags columns)")
    parser.add_argument("--energy-threshold", type=float, default=8.0)
    parser.add_argument("--ax", type=float, default=0.25)
    parser.add_argument("--window-states", type=int, default=10)
    parser.add_argument("--std2-bin", default=os.path.expanduser("~/local_bin/std2"))
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    threshold, matched, mismatched = derive_threshold(args.calibration_csv)
    print(f"Calibration: {len(matched)} character-matched conformers "
          f"(contrast range {min(matched):.1f}-{max(matched):.1f}x), "
          f"{len(mismatched)} character-mismatched "
          f"(contrast range {min(mismatched):.1f}-{max(mismatched):.1f}x)"
          if mismatched else
          f"Calibration: {len(matched)} character-matched conformers, "
          f"0 mismatched -- no threshold can be derived (nothing wrong "
          f"was ever observed in the labeled set).")
    if threshold is not None:
        n_matched_below = sum(1 for c in matched if c <= threshold)
        print(f"Threshold = {threshold:.2f}x (highest contrast among known "
              f"mismatches). Note: {n_matched_below}/{len(matched)} known-"
              f"CORRECT picks also fall at or below this threshold -- "
              f"the two distributions overlap, so this flag is a "
              f"conservative screen, not a clean separator.")

    with open(args.job_list) as f:
        job_dirs = [line.strip() for line in f if line.strip()]

    rows = []
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
        if states is None:
            print(f"WARNING: could not parse std2 output for {jd}")
            continue
        window = sorted(states, key=lambda s: s["energy_eV"])[:args.window_states]
        max_f_state = max(window, key=lambda s: s["osc"])
        other_f = sorted((s["osc"] for s in window if s is not max_f_state), reverse=True)
        second_f = other_f[0] if other_f else 0.0
        contrast = (max_f_state["osc"] / second_f) if second_f > 1e-9 else float("inf")
        is_uncertain = (threshold is not None and contrast != float("inf")
                         and contrast <= threshold)

        rows.append({
            "job_dir": jd, "s1_energy_eV_stda": max_f_state["energy_eV"],
            "s1_osc_stda": max_f_state["osc"],
            "f_contrast_ratio": contrast if contrast != float("inf") else "inf",
            "uncertain_selection": is_uncertain,
        })
        print(f"{jd}: E={max_f_state['energy_eV']:.4f} eV, "
              f"f={max_f_state['osc']:.5f}, contrast="
              f"{contrast:.2f}x" if contrast != float("inf") else "contrast=inf")

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["job_dir", "s1_energy_eV_stda", "s1_osc_stda",
                           "f_contrast_ratio", "uncertain_selection"])
        writer.writeheader()
        writer.writerows(rows)

    n_uncertain = sum(1 for r in rows if r["uncertain_selection"])
    print(f"\n=== Summary (N={len(rows)}) ===")
    print(f"Flagged 'uncertain_selection': {n_uncertain}/{len(rows)} "
          f"({100*n_uncertain/len(rows):.1f}%)" if rows else "No conformers processed.")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
