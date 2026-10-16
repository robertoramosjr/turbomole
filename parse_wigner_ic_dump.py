#!/usr/bin/env python3
"""
parse_wigner_ic_dump.py

Consolidates a Newton-X initial-condition batch (one folder per split job,
e.g. I1..I100 from split_initcond.pl) into a single flat dataset covering
every (Wigner point, excited state) pair -- not just the states Newton-X
itself extracted into initial_condition.<nis>.<nfs>.

Each batch folder is expected to contain:
  - TEMP/exspectrum: one "# Excitation spectrum ..." block per point
    (equilibrium geometry first, then each Wigner-sampled point in order),
    each block listing all computed states (energy, oscillator strengths).
  - initcond.log: one [STATE INFORMATION]...[END STATE INFORMATION]
    section per point, in the same order, each containing one "Dominant
    contributions:" block per state (also in state order) with the
    leading (highest-weight) occ->virt orbital pair for that state.

The two sources are joined by POSITION (Nth exspectrum block <-> Nth
[STATE INFORMATION] section; Nth state row <-> Nth "Dominant
contributions" block within that section), since neither file prints an
explicit Wigner-point index. This assumes both files were written by the
same run without exspectrum being reused/appended-to across a different
run of the same folder (see --check-consistency for a sanity check
against the folder's own initial_condition.<nis>.<nfs>, if present).

Output: one row per (point, state), plus a separate equilibrium-geometry
cross-check dataset (one row per batch folder, since every folder
independently recomputes the same equilibrium point).

Usage:
    python ~/scripts/parse_wigner_ic_dump.py \
        --batch-dir INITIAL_CONDITIONS --pattern "I*" \
        --states-per-point 10 \
        --output wigner_1000_states_dump.csv \
        --equilibrium-output wigner_equilibrium_crosscheck.csv
"""

import argparse
import glob
import os
import re

EXSPECTRUM_HEADER_RE = re.compile(r"^#\s*Excitation spectrum")
EXSPECTRUM_ROW_RE = re.compile(
    r"^\s*(\d+)\s+a\s+([-\d.]+)\s+([-\d.]+)\s+[-\d.D+]+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$"
)

DOMINANT_HEADER_RE = re.compile(r"^\s*Dominant contributions:\s*$")
DOMINANT_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(\S+)\s+([-\d.]+)\s+(\d+)\s+(\S+)\s+([-\d.]+)\s+([\d.]+)\s*$"
)
STATE_INFO_START_RE = re.compile(r"^\[STATE INFORMATION\]\s*$")
STATE_INFO_END_RE = re.compile(r"^\[END STATE INFORMATION\]\s*$")


def parse_exspectrum(path, states_per_point):
    """Returns a list of point-blocks; each is a list of
    `states_per_point` dicts with energy/oscillator-strength data."""
    points = []
    current = []
    with open(path, "r") as f:
        for line in f:
            if EXSPECTRUM_HEADER_RE.match(line):
                if current:
                    points.append(current)
                current = []
                continue
            match = EXSPECTRUM_ROW_RE.match(line)
            if match:
                state, e_eh, e_ev, e_nm, osc_vel, osc_len = match.groups()
                current.append({
                    "state": int(state), "energy_eV": float(e_ev),
                    "energy_nm": float(e_nm), "osc_vel": float(osc_vel),
                    "osc_len": float(osc_len),
                })
    if current:
        points.append(current)

    bad = [i for i, p in enumerate(points, 1) if len(p) != states_per_point]
    if bad:
        raise ValueError(
            f"{path}: point block(s) {bad} have != {states_per_point} state rows "
            f"-- check --states-per-point or an incomplete/truncated exspectrum.")
    return points


def parse_dominant_contributions(path, states_per_point):
    """Returns a list of point-blocks (same order as exspectrum); each is
    a list of `states_per_point` dicts with the leading occ->virt pair."""
    points = []
    current = []
    in_state_info = False
    awaiting_header = False
    awaiting_row = False
    with open(path, "r") as f:
        for line in f:
            if STATE_INFO_START_RE.match(line):
                in_state_info = True
                current = []
                continue
            if STATE_INFO_END_RE.match(line):
                in_state_info = False
                points.append(current)
                continue
            if not in_state_info:
                continue

            if DOMINANT_HEADER_RE.match(line):
                awaiting_header = True
                continue
            if awaiting_header:
                if line.strip():
                    awaiting_header = False
                    awaiting_row = True
                continue
            if awaiting_row:
                if not line.strip():
                    continue
                match = DOMINANT_ROW_RE.match(line)
                awaiting_row = False
                if match:
                    occ_idx, occ_sym, occ_e, virt_idx, virt_sym, virt_e, weight = match.groups()
                    current.append({
                        "orbital": f"{occ_idx}{occ_sym}->{virt_idx}{virt_sym}",
                        "weight_pct": float(weight),
                    })

    bad = [i for i, p in enumerate(points, 1) if len(p) != states_per_point]
    if bad:
        raise ValueError(
            f"{path}: [STATE INFORMATION] block(s) {bad} have != {states_per_point} "
            f"'Dominant contributions' entries -- check --states-per-point.")
    return points


def main():
    parser = argparse.ArgumentParser(
        description="Consolidate a split Newton-X initial-condition batch (exspectrum + "
                    "initcond.log per folder) into one flat per-(point,state) dataset.")
    parser.add_argument("--batch-dir", required=True,
                         help="Directory containing the split batch folders (e.g. INITIAL_CONDITIONS)")
    parser.add_argument("--pattern", default="I[0-9]*",
                         help="Glob pattern (within --batch-dir) for the batch folders")
    parser.add_argument("--exspectrum-rel", default="TEMP/exspectrum",
                         help="Path to exspectrum relative to each batch folder")
    parser.add_argument("--initcond-log", default="initcond.log",
                         help="Name of the Newton-X log file inside each batch folder")
    parser.add_argument("--states-per-point", type=int, required=True,
                         help="Number of excited states requested per point (Turbomole $soes)")
    parser.add_argument("--output", default="wigner_states_dump.csv")
    parser.add_argument("--equilibrium-output", default="wigner_equilibrium_crosscheck.csv")
    args = parser.parse_args()

    folders = sorted(
        glob.glob(os.path.join(args.batch_dir, args.pattern)),
        key=lambda p: int(re.search(r"\d+", os.path.basename(p)).group()),
    )
    if not folders:
        raise ValueError(f"No folders matched {args.pattern} under {args.batch_dir}")

    state_rows = []
    equilibrium_rows = []
    point_id = 0

    for folder in folders:
        folder_name = os.path.basename(folder)
        exspectrum_path = os.path.join(folder, args.exspectrum_rel)
        initcond_log_path = os.path.join(folder, args.initcond_log)

        energy_points = parse_exspectrum(exspectrum_path, args.states_per_point)
        orbital_points = parse_dominant_contributions(initcond_log_path, args.states_per_point)

        if len(energy_points) != len(orbital_points):
            raise ValueError(
                f"{folder_name}: {len(energy_points)} exspectrum point-blocks but "
                f"{len(orbital_points)} [STATE INFORMATION] sections -- can't align by position.")

        # First block in each folder is always the equilibrium geometry
        # (0 velocity), recomputed independently by every split folder.
        equilibrium_energy, equilibrium_orbitals = energy_points[0], orbital_points[0]
        for state_energy, state_orbital in zip(equilibrium_energy, equilibrium_orbitals):
            equilibrium_rows.append({
                "folder": folder_name, "state": state_energy["state"],
                "energy_eV": state_energy["energy_eV"], "energy_nm": state_energy["energy_nm"],
                "osc_len": state_energy["osc_len"], "osc_vel": state_energy["osc_vel"],
                "dominant_orbital": state_orbital["orbital"],
                "weight_pct": state_orbital["weight_pct"],
            })

        for wigner_index, (energy_block, orbital_block) in enumerate(
                zip(energy_points[1:], orbital_points[1:]), start=1):
            point_id += 1
            for state_energy, state_orbital in zip(energy_block, orbital_block):
                state_rows.append({
                    "point_id": point_id, "wigner_index": wigner_index, "folder": folder_name,
                    "state": state_energy["state"], "energy_eV": state_energy["energy_eV"],
                    "energy_nm": state_energy["energy_nm"], "osc_len": state_energy["osc_len"],
                    "osc_vel": state_energy["osc_vel"],
                    "dominant_orbital": state_orbital["orbital"],
                    "weight_pct": state_orbital["weight_pct"],
                })

    print(f"Parsed {len(folders)} batch folders, {point_id} Wigner points, "
          f"{args.states_per_point} states each -> {len(state_rows)} rows "
          f"(expected {point_id * args.states_per_point})")

    fieldnames = ["point_id", "wigner_index", "folder", "state", "energy_eV", "energy_nm",
                  "osc_len", "osc_vel", "dominant_orbital", "weight_pct"]
    with open(args.output, "w") as f:
        f.write(",".join(fieldnames) + "\n")
        for r in state_rows:
            f.write(",".join(
                f"{r[k]:.6f}" if isinstance(r[k], float) else str(r[k])
                for k in fieldnames) + "\n")
    print(f"Per-(point,state) dataset saved -> {args.output}")

    eq_fieldnames = ["folder", "state", "energy_eV", "energy_nm", "osc_len", "osc_vel",
                      "dominant_orbital", "weight_pct"]
    with open(args.equilibrium_output, "w") as f:
        f.write(",".join(eq_fieldnames) + "\n")
        for r in equilibrium_rows:
            f.write(",".join(
                f"{r[k]:.6f}" if isinstance(r[k], float) else str(r[k])
                for k in eq_fieldnames) + "\n")
    print(f"Equilibrium-geometry cross-check dataset saved -> {args.equilibrium_output} "
          f"({len(folders)} folders x {args.states_per_point} states)")


if __name__ == "__main__":
    main()
