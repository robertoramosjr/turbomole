#!/usr/bin/env python3
"""
screen_confident_conformers_for_hidden_state_risk.py

Screens the neutral state's "confident" (f-contrast > threshold)
conformers for the NEW failure signature confirmed by job_0023's DFT
confirmation (Sec. 9.23-9.24 of protonation_and_conformer_investigation.md):
sTDA can be internally CONFIDENT (high contrast) in picking the
canonical HOMO->LUMO (offset (0,0)) state as the brightest in its
window, while HSE06's real S1 is a different, non-canonical state that
sTDA ranks far lower and dimmer -- a failure mode the f-contrast check
cannot see at all, since it only measures how much the CHOSEN state
dominates its own window, not whether that state corresponds to
HSE06's real answer.

For job_0023, the tell in hindsight was that a non-canonical
(offset != (0,0)) state was still PRESENT somewhere in the low-energy
sTDA window with non-negligible oscillator strength (f=0.026 at rank
5, well above the numerical noise floor of the dimmest window states,
~0.0005-0.001) even though it was nowhere near competitive with the
chosen state's f=0.626. This script's risk signature is exactly that:
does ANY non-canonical state exist in the window with f at or above
--risk-f-threshold, regardless of how far behind the chosen state it
is. This is a broad, low-specificity screen (a real molecule's low-
lying spectrum can have weak non-canonical contamination without that
implying a real problem) -- it flags candidates for a closer look or
for inclusion in a future confirmation batch, not confirmed errors.

Usage:
    python ~/work_turbomole/scripts/screen_confident_conformers_for_hidden_state_risk.py \
        --job-list confident_job_dirs.txt \
        --r2scan3c-jobs layer_a/conformer_jobs \
        --ax 0.25 --window-states 10 --risk-f-threshold 0.02 \
        --output layer_b/hidden_state_risk_screen.csv
"""

import argparse
import os
import re
import subprocess
import tempfile
import csv

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
    return occ_l - local_homo, virt_l - (local_homo + 1)


def main():
    parser = argparse.ArgumentParser(
        description="Screen conformers for a non-canonical low-lying sTDA "
                    "state with non-negligible f, present regardless of the "
                    "chosen (max-f) state's own contrast -- the signature "
                    "behind job_0023's confirmed 16x oscillator-strength miss.")
    parser.add_argument("--job-list", required=True)
    parser.add_argument("--r2scan3c-jobs", required=True)
    parser.add_argument("--energy-threshold", type=float, default=8.0)
    parser.add_argument("--ax", type=float, default=0.25)
    parser.add_argument("--window-states", type=int, default=10)
    parser.add_argument("--risk-f-threshold", type=float, default=0.02,
                         help="A non-canonical state at or above this f is "
                              "flagged as a hidden-state risk [default: "
                              "0.02, matching job_0023's flagged state, "
                              "f=0.026, comfortably above the ~0.0005-0.001 "
                              "noise floor of the dimmest window states]")
    parser.add_argument("--std2-bin", default=os.path.expanduser("~/local_bin/std2"))
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

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
        local_homo = parse_local_homo_index(stdout)
        if states is None or local_homo is None:
            print(f"WARNING: could not parse std2 output for {jd}")
            continue

        window = sorted(states, key=lambda s: s["energy_eV"])[:args.window_states]
        chosen = max(window, key=lambda s: s["osc"])
        chosen_offset = dominant_offset(chosen, local_homo)

        risky_states = []
        for s in window:
            if s is chosen:
                continue
            off = dominant_offset(s, local_homo)
            if off != (0, 0) and s["osc"] >= args.risk_f_threshold:
                risky_states.append((off, s["osc"], s["energy_eV"]))

        flagged = len(risky_states) > 0
        risky_str = "; ".join(f"offset={o} f={f:.4f} E={e:.3f}eV"
                                for o, f, e in risky_states)
        rows.append({
            "job_dir": jd, "chosen_offset": str(chosen_offset),
            "chosen_f": chosen["osc"], "hidden_state_risk": flagged,
            "risky_states": risky_str,
        })
        print(f"{jd}: chosen_offset={chosen_offset}, f={chosen['osc']:.4f} -- "
              f"{'HIDDEN_STATE_RISK: ' + risky_str if flagged else 'clean'}")

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["job_dir", "chosen_offset", "chosen_f",
                           "hidden_state_risk", "risky_states"])
        writer.writeheader()
        writer.writerows(rows)

    n_flagged = sum(1 for r in rows if r["hidden_state_risk"])
    print(f"\n=== Summary (N={len(rows)}) ===")
    print(f"Flagged hidden_state_risk: {n_flagged}/{len(rows)} "
          f"({100 * n_flagged / len(rows):.1f}%)" if rows else "No conformers processed.")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
