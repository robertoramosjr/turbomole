#!/usr/bin/env python3
"""
analyze_stda_state_selection.py

Validates the "select the sTDA state with the largest oscillator
strength within the energy window" criterion (proposed to replace
"take state 1 by energy rank", after the ax-grid diagnostic showed
that sTDA robustly places a spurious dark state below the true bright
S1-like transition for the neutral and monoanion_b protonation states,
across the ENTIRE ax=0.0-0.30 range -- see
protonation_and_conformer_investigation.md, oscillator-strength
diagnostic sections).

For every conformer in --job-list, runs tm2molden + std2 once (same
pipeline as batch_stda_s1.py) and reports three checks together:

  1. Energy RANK of the state with the largest oscillator strength
     among all states std2 prints (not assumed to be state 1).
  2. CONTRAST: ratio of that state's f to the second-largest f among
     the other states -- a low ratio means "largest f" is an ambiguous
     criterion for that conformer.
  3. ORBITAL CHARACTER: whether the max-f state's dominant occupied-
     >virtual contribution is the local HOMO->LUMO transition, in the
     SAME local numbering std2 itself reports for that run (std2
     truncates to an active MO window and renumbers locally per run,
     so absolute indices are not comparable across ax values or to
     Turbomole's own numbering -- "oMOs in TDA: N" in std2's stdout
     gives the local HOMO index N directly). Compared against the
     reference HSE06 escf.out's own "Dominant contributions" table for
     S1, using Turbomole's own "number of occupied orbitals" line to
     identify ITS HOMO index independently.

Usage:
    python ~/work_turbomole/scripts/analyze_stda_state_selection.py \
        --job-list job_dirs.txt \
        --r2scan3c-jobs layer_a/conformer_jobs \
        --labeled-jobs layer_b/labeled_jobs \
        --ax 0.25 \
        --output layer_b/stda_state_selection_check.csv
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
    s1_energy = s1_osc = dom_occ = dom_virt = dom_weight = None
    for i, line in enumerate(lines):
        if "1 singlet a excitation" in line:
            for j in range(i, min(i + 60, len(lines))):
                if "Excitation energy / eV:" in lines[j]:
                    s1_energy = float(lines[j].split(":")[1].strip())
                if "length representation:" in lines[j] and s1_osc is None \
                        and "Oscillator" in "\n".join(lines[i:j]):
                    try:
                        s1_osc = float(lines[j].split(":")[1].strip())
                    except ValueError:
                        pass
                if "Dominant contributions" in lines[j]:
                    for k in range(j, min(j + 6, len(lines))):
                        dm = HSE06_DOMINANT_RE.match(lines[k])
                        if dm:
                            dom_occ, dom_virt = int(dm.group(1)), int(dm.group(2))
                            dom_weight = float(dm.group(3))
                            break
                    break
            break
    return {
        "homo_index": homo_index, "s1_energy_eV": s1_energy, "s1_osc": s1_osc,
        "dom_occ": dom_occ, "dom_virt": dom_virt, "dom_weight_pct": dom_weight,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Check rank, f-contrast, and orbital character of the "
                    "sTDA state selected by largest oscillator strength, "
                    "against the HSE06 reference S1.")
    parser.add_argument("--job-list", required=True)
    parser.add_argument("--r2scan3c-jobs", required=True)
    parser.add_argument("--labeled-jobs", required=True,
                         help="Folder of job_NNNN/escf.out (HSE06 production labels)")
    parser.add_argument("--energy-threshold", type=float, default=8.0)
    parser.add_argument("--ax", type=float, default=0.25,
                         help="[default: 0.25 -- matches HSE06's exact-exchange "
                              "fraction; the ax-grid diagnostic showed this gives "
                              "a cleaner, less-mixed bright-state character than "
                              "ax=0.0, even though it does not fix state ranking]")
    parser.add_argument("--std2-bin", default=os.path.expanduser("~/local_bin/std2"))
    parser.add_argument("--contrast-flag-threshold", type=float, default=10.0)
    parser.add_argument("--rank-flag-threshold", type=int, default=5)
    parser.add_argument("--max-window-states", type=int, default=None,
                         help="If set, restrict the max-f search (and the "
                              "contrast comparison) to only the first N "
                              "states by energy, instead of the full sTDA "
                              "spectral window -- guards against a spurious "
                              "high-energy state (large f, wrong orbital "
                              "character) hijacking the selection when the "
                              "true bright state sits at a low rank anyway.")
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

        hse06 = parse_hse06_s1(escf_path)

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

        states_sorted_by_energy = sorted(states, key=lambda s: s["energy_eV"])
        window = states_sorted_by_energy[:args.max_window_states] \
            if args.max_window_states else states_sorted_by_energy
        max_f_state = max(window, key=lambda s: s["osc"])
        rank = states_sorted_by_energy.index(max_f_state) + 1

        other_f = sorted(
            (s["osc"] for s in window if s is not max_f_state),
            reverse=True,
        )
        second_f = other_f[0] if other_f else 0.0
        contrast = (max_f_state["osc"] / second_f) if second_f > 1e-9 else float("inf")

        dom = max(max_f_state["contribs"], key=lambda c: abs(c[0])) \
            if max_f_state["contribs"] else (None, None, None)
        dom_coeff, dom_occ_local, dom_virt_local = dom
        is_homo_lumo = (dom_occ_local == local_homo and dom_virt_local == local_homo + 1)

        flags = []
        if rank > args.rank_flag_threshold:
            flags.append("RANK_OUT_OF_WINDOW")
        if contrast < args.contrast_flag_threshold:
            flags.append("LOW_CONTRAST")
        if not is_homo_lumo:
            flags.append("CHARACTER_MISMATCH")

        row = {
            "job_dir": jd,
            "hse06_s1_energy_eV": hse06["s1_energy_eV"],
            "hse06_s1_osc": hse06["s1_osc"],
            "hse06_homo_index": hse06["homo_index"],
            "hse06_dom_occ": hse06["dom_occ"],
            "hse06_dom_virt": hse06["dom_virt"],
            "hse06_dom_weight_pct": hse06["dom_weight_pct"],
            "hse06_is_homo_lumo": hse06["dom_occ"] == hse06["homo_index"]
                and hse06["dom_virt"] == hse06["homo_index"] + 1,
            "n_states_parsed": len(states),
            "stda_max_f_rank": rank,
            "stda_max_f_energy_eV": max_f_state["energy_eV"],
            "stda_max_f_osc": max_f_state["osc"],
            "stda_f_contrast_ratio": contrast,
            "stda_local_homo_index": local_homo,
            "stda_dom_occ_local": dom_occ_local,
            "stda_dom_virt_local": dom_virt_local,
            "stda_dom_coeff": dom_coeff,
            "stda_is_homo_lumo_match": is_homo_lumo,
            "flags": ";".join(flags) if flags else "",
        }
        rows.append(row)
        print(f"{jd}: rank={rank}, contrast={contrast:.1f}x, "
              f"HOMO->LUMO match={is_homo_lumo}"
              + (f"  [{row['flags']}]" if flags else ""))

    if not rows:
        print("No conformers processed.")
        return

    fieldnames = list(rows[0].keys())
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    ranks = [r["stda_max_f_rank"] for r in rows]
    contrasts = [r["stda_f_contrast_ratio"] for r in rows if r["stda_f_contrast_ratio"] != float("inf")]
    n_rank_flag = sum(1 for r in rows if "RANK_OUT_OF_WINDOW" in r["flags"])
    n_contrast_flag = sum(1 for r in rows if "LOW_CONTRAST" in r["flags"])
    n_char_flag = sum(1 for r in rows if "CHARACTER_MISMATCH" in r["flags"])
    n_hse06_not_homo_lumo = sum(1 for r in rows if not r["hse06_is_homo_lumo"])

    print(f"\n=== Summary (N={len(rows)}) ===")
    print(f"Rank distribution: min={min(ranks)}, max={max(ranks)}, "
          f"rank>{args.rank_flag_threshold} in {n_rank_flag}/{len(rows)}")
    if contrasts:
        print(f"Contrast ratio: min={min(contrasts):.1f}x, "
              f"median={sorted(contrasts)[len(contrasts)//2]:.1f}x, "
              f"<{args.contrast_flag_threshold}x in {n_contrast_flag}/{len(rows)}")
    print(f"HOMO->LUMO character mismatch (sTDA vs. own HSE06 reference): "
          f"{n_char_flag}/{len(rows)}")
    if n_hse06_not_homo_lumo:
        print(f"NOTE: {n_hse06_not_homo_lumo}/{len(rows)} conformers have an "
              f"HSE06 S1 that is NOT itself dominant HOMO->LUMO -- character "
              f"check for those relies on the printed dominant pair, not the "
              f"HOMO->LUMO shortcut.")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
