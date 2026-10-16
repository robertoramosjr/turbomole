#!/usr/bin/env python3
"""
parse_s1_scan.py

Extracts the results of a relaxed excited-state (egrad) dihedral scan
built with Turbomole's `scanprep -c` (concatenated scan: each point
starts from the previous point's already-relaxed geometry). Reads each
tors_<angle> folder and reports, for the tracked state ($exopt) and its
neighbor:
  - the actual CNNC dihedral in the final geometry (sanity check against
    the requested scan value)
  - E(S0), the tracked state's excitation energy (gap) and oscillator
    strength, from the last-written $exspectrum block
  - the dominant occ->virt orbital pair (and its weight) for both the
    tracked state and the next one, from job.last's "Dominant
    contributions" blocks -- this is what lets you tell whether a state
    is still n->pi* or has swapped character to pi->pi* along the scan,
    instead of trusting the S1/S2 index blindly.

Points with no GEO_OPT_CONVERGED marker are reported as failed, with
whatever partial information survives (last successful gradient cycle's
excitation energy, if any) rather than silently dropped.

Usage:
    python ~/scripts/parse_s1_scan.py --scan-dir scan_s1_b3lyp \
        --dihedral-atoms 5,23,24,8 \
        --tracked-state 1 --neighbor-state 2 \
        --output s1_scan_b3lyp.dat
"""

import argparse
import glob
import os
import re

import numpy as np

EXSPECTRUM_HEADER_RE = re.compile(r"^#\s*Excitation spectrum")
EXSPECTRUM_ROW_RE = re.compile(
    r"^\s*(\d+)\s+a\s+([-\d.]+)\s+([-\d.]+)\s+[-\d.D+]+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$"
)
STATE_HEADER_RE = re.compile(r"^\s*(\d+)\s+singlet\s+\S+\s+excitation\s*$")
DOMINANT_HEADER_RE = re.compile(r"^\s*Dominant contributions:\s*$")
DOMINANT_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(\S+)\s+([-\d.]+)\s+(\d+)\s+(\S+)\s+([-\d.]+)\s+([\d.]+)\s*$"
)


def read_coord_dihedral(coord_path, atom_indices):
    atoms = []
    with open(coord_path, "r") as f:
        reading = False
        for line in f:
            if line.startswith("$coord"):
                reading = True
                continue
            if reading and line.startswith("$"):
                break
            if reading:
                parts = line.split()
                atoms.append([float(parts[0]), float(parts[1]), float(parts[2])])
    pos = np.array(atoms)
    i0, i1, i2, i3 = [i - 1 for i in atom_indices]
    p0, p1, p2, p3 = pos[i0], pos[i1], pos[i2], pos[i3]
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def last_exspectrum_block(path):
    blocks = []
    current = None
    with open(path, "r") as f:
        for line in f:
            if EXSPECTRUM_HEADER_RE.match(line):
                if current:
                    blocks.append(current)
                current = {}
                continue
            match = EXSPECTRUM_ROW_RE.match(line)
            if match and current is not None:
                state, e_eh, e_ev, e_nm, osc_vel, osc_len = match.groups()
                current[int(state)] = {
                    "energy_eV": float(e_ev), "energy_nm": float(e_nm),
                    "osc_vel": float(osc_vel), "osc_len": float(osc_len),
                }
    if current:
        blocks.append(current)
    return blocks[-1] if blocks else None


def last_dominant_contributions(path):
    """Returns {state: {occ_orbital, virt_orbital, weight_pct}} using the
    LAST occurrence of each state's block in job.last (i.e. the final,
    converged-geometry values, since job.last accumulates every SCF/egrad
    call made during the optimization)."""
    contributions = {}
    current_state = None
    awaiting_header = False
    awaiting_row = False
    with open(path, "r") as f:
        for line in f:
            header_match = STATE_HEADER_RE.match(line)
            if header_match:
                current_state = int(header_match.group(1))
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
                if match and current_state is not None:
                    occ_idx, occ_sym, _occ_e, virt_idx, virt_sym, _virt_e, weight = match.groups()
                    contributions[current_state] = {
                        "orbital": f"{occ_idx}{occ_sym}->{virt_idx}{virt_sym}",
                        "weight_pct": float(weight),
                    }
    return contributions


def last_ground_state_energy(energy_path):
    last_scf = None
    with open(energy_path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                last_scf = float(parts[1])
    return last_scf


def angle_from_folder(folder):
    return float(os.path.basename(folder.rstrip("/")).replace("tors_", ""))


def main():
    parser = argparse.ArgumentParser(
        description="Extract energetics and orbital character from a relaxed "
                    "excited-state dihedral scan (Turbomole scanprep -c + egrad).")
    parser.add_argument("--scan-dir", required=True,
                         help="Directory containing tors_<angle> subfolders")
    parser.add_argument("--dihedral-atoms", required=True,
                         help="Comma-separated 1-based atom indices for the scanned "
                              "dihedral, e.g. 5,23,24,8")
    parser.add_argument("--tracked-state", type=int, required=True,
                         help="State index tracked by $exopt (e.g. 1 for S1)")
    parser.add_argument("--neighbor-state", type=int, required=True,
                         help="Adjacent state to also report for character-swap checks")
    parser.add_argument("--output", default="s1_scan.dat")
    args = parser.parse_args()

    dihedral_atoms = [int(x) for x in args.dihedral_atoms.split(",")]
    folders = sorted(
        (p for p in glob.glob(os.path.join(args.scan_dir, "tors_*")) if os.path.isdir(p)),
        key=angle_from_folder, reverse=True)

    rows = []
    for folder in folders:
        target_angle = angle_from_folder(folder)
        converged = os.path.exists(os.path.join(folder, "GEO_OPT_CONVERGED"))
        exspectrum_path = os.path.join(folder, "exspectrum")
        coord_path = os.path.join(folder, "coord")
        job_last_path = os.path.join(folder, "job.last")
        energy_path = os.path.join(folder, "energy")

        row = {"target_cnnc": target_angle, "converged": converged,
               "actual_cnnc": "", "e_s0_Eh": "", "gap_eV": "", "gap_nm": "",
               "osc_len": "", "tracked_orbital": "", "tracked_weight_pct": "",
               "neighbor_gap_eV": "", "neighbor_orbital": "", "neighbor_weight_pct": "",
               "note": ""}

        if os.path.exists(coord_path):
            try:
                row["actual_cnnc"] = read_coord_dihedral(coord_path, dihedral_atoms)
            except Exception as e:
                row["note"] += f"coord parse failed: {e}; "

        exspectrum = last_exspectrum_block(exspectrum_path) if os.path.exists(exspectrum_path) else None
        if exspectrum and args.tracked_state in exspectrum:
            row["gap_eV"] = exspectrum[args.tracked_state]["energy_eV"]
            row["gap_nm"] = exspectrum[args.tracked_state]["energy_nm"]
            row["osc_len"] = exspectrum[args.tracked_state]["osc_len"]
        if exspectrum and args.neighbor_state in exspectrum:
            row["neighbor_gap_eV"] = exspectrum[args.neighbor_state]["energy_eV"]

        if os.path.exists(energy_path):
            e_s0 = last_ground_state_energy(energy_path)
            if e_s0 is not None:
                row["e_s0_Eh"] = e_s0

        if os.path.exists(job_last_path):
            contributions = last_dominant_contributions(job_last_path)
            if args.tracked_state in contributions:
                row["tracked_orbital"] = contributions[args.tracked_state]["orbital"]
                row["tracked_weight_pct"] = contributions[args.tracked_state]["weight_pct"]
            if args.neighbor_state in contributions:
                row["neighbor_orbital"] = contributions[args.neighbor_state]["orbital"]
                row["neighbor_weight_pct"] = contributions[args.neighbor_state]["weight_pct"]

        if not converged:
            failed_marker = os.path.join(folder, "GEO_OPT_FAILED")
            if os.path.exists(failed_marker):
                with open(failed_marker) as f:
                    row["note"] += "FAILED: " + f.readline().strip() + "; "
            else:
                row["note"] += "no GEO_OPT_CONVERGED/GEO_OPT_FAILED marker found; "
            if row["gap_eV"] == "":
                row["note"] += "no exspectrum data recovered; "

        rows.append(row)

    fieldnames = ["target_cnnc", "actual_cnnc", "converged", "e_s0_Eh", "gap_eV", "gap_nm",
                  "osc_len", "tracked_orbital", "tracked_weight_pct",
                  "neighbor_gap_eV", "neighbor_orbital", "neighbor_weight_pct", "note"]
    with open(args.output, "w") as f:
        f.write("descriptor " + ",".join(fieldnames) + "\n")
        for r in rows:
            values = []
            for k in fieldnames:
                v = r[k]
                if isinstance(v, float):
                    v = f"{v:.6f}"
                elif isinstance(v, bool):
                    v = "1" if v else "0"
                elif v == "":
                    v = "NA"
                values.append(str(v).replace(" ", "_"))
            f.write(" ".join(values) + "\n")

    n_converged = sum(1 for r in rows if r["converged"])
    print(f"{len(rows)} scan points, {n_converged} converged, {len(rows)-n_converged} failed")
    for r in rows:
        status = "OK" if r["converged"] else "FAILED"
        print(f"  CNNC={r['target_cnnc']:.0f} [{status}]: "
              f"actual_cnnc={r['actual_cnnc']}, gap={r['gap_eV']} eV, "
              f"f={r['osc_len']}, tracked={r['tracked_orbital']} ({r['tracked_weight_pct']}%)"
              f"{'  ' + r['note'] if r['note'] else ''}")
    print(f"\nDataset saved -> {args.output}")


if __name__ == "__main__":
    main()
