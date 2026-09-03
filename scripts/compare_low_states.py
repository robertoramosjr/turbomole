#!/usr/bin/env python3
"""
compare_low_states.py

Extracts the N lowest excited states from multiple escf.out files --
e.g. TD-DFT runs with different functionals and a GW-BSE/TDA run --
and builds a single long-format comparison table, state by state,
across all runs.

For each state, two different energies are reported:
  - energy_eV: the actual (many-body) excitation energy from the
    SUMMARY OF EXCITATION ENERGIES table -- includes electron-hole
    interaction / kernel response, NOT just an orbital difference.
  - orbital_gap_eV: the single-particle energy difference between the
    dominant occupied and virtual orbital of that state's "Dominant
    contributions" block (i.e. occ_energy_eV / virt_energy_eV). This
    coincides with the true HOMO-LUMO gap only for states whose
    dominant transition actually is HOMO->LUMO -- for other states it
    is the gap of whichever occ/virt orbital pair dominates that
    state, so occ_orbital/virt_orbital are reported alongside it.
  - binding_energy_eV: orbital_gap_eV - energy_eV, i.e. the exciton
    binding energy computed against the DOMINANT transition's own
    orbital pair (not necessarily HOMO-LUMO). This is a different,
    per-state quantity from compute_exciton_binding.py, which always
    uses the fixed HOMO/LUMO pair -- do not mix the two up.

The excitation table format (SUMMARY OF EXCITATION ENERGIES AND
DIPOLE OSCILLATOR STRENGTHS, and the per-state "N singlet a
excitation" / "Dominant contributions" blocks) is identical between
TD-DFT/RPA and BSE/TDA escf runs, so the same regexes cover both.

Usage:
    python ~/scripts/compare_low_states.py \
        --run hse06-tddft:artepillin_C_d3_hse06/escf_d3.out \
        --run b3lyp-tddft:artepillin_C_d3_b3lyp/escf_b3lyp.out \
        --run camb3lyp-tddft:artepillin_C_d3_camb3lyp/escf_camb3lyp.out \
        --run hse06-bse-tda:artepillin_C_d3_hse06_gwbse/escf_bse.out \
        --nstates 10 \
        --output low_states_comparison.csv
"""

import argparse
import re

TABLE_ROW_RE = re.compile(
    r"^\s*\|\s*(\d+)\s+a\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
)

STATE_HEADER_RE = re.compile(r"^\s*(\d+)\s+singlet\s+\S+\s+excitation\s*$")
DOMINANT_HEADER_RE = re.compile(r"^\s*Dominant contributions:\s*$")
DOMINANT_ROW_RE = re.compile(
    r"^\s*(\d+)\s+(\S+)\s+([-\d.]+)\s+(\d+)\s+(\S+)\s+([-\d.]+)\s+([\d.]+)\s*$"
)


def parse_escf_table(escf_path):
    rows = []
    with open(escf_path, "r") as f:
        for line in f:
            match = TABLE_ROW_RE.match(line)
            if match:
                state, e_eh, e_ev, e_cm1, e_nm, osc_vel, osc_len = match.groups()
                rows.append({
                    "state": int(state),
                    "energy_eV": float(e_ev),
                    "energy_nm": float(e_nm),
                    "osc_vel": float(osc_vel),
                    "osc_len": float(osc_len),
                })
    if not rows:
        raise ValueError(f"No excitation rows found in {escf_path}")
    return sorted(rows, key=lambda r: r["state"])


def parse_dominant_transitions(escf_path):
    """Dominant (highest-weight) occ->virt orbital pair per state.

    Only the first data row under each "Dominant contributions:"
    block is kept -- the rows are printed in descending |coeff|^2
    order, so the first one is the dominant single-particle
    transition for that state.
    """
    transitions = {}
    current_state = None
    awaiting_column_header = False
    awaiting_data_row = False
    with open(escf_path, "r") as f:
        for line in f:
            header_match = STATE_HEADER_RE.match(line)
            if header_match:
                current_state = int(header_match.group(1))
                continue
            if DOMINANT_HEADER_RE.match(line):
                awaiting_column_header = True
                continue
            if awaiting_column_header:
                if line.strip():
                    awaiting_column_header = False
                    awaiting_data_row = True
                continue
            if awaiting_data_row:
                if not line.strip():
                    continue
                row_match = DOMINANT_ROW_RE.match(line)
                awaiting_data_row = False
                if row_match and current_state is not None:
                    occ_idx, occ_sym, occ_e, virt_idx, virt_sym, virt_e, weight = row_match.groups()
                    transitions[current_state] = {
                        "occ_orbital": f"{occ_idx}{occ_sym}",
                        "occ_energy_eV": float(occ_e),
                        "virt_orbital": f"{virt_idx}{virt_sym}",
                        "virt_energy_eV": float(virt_e),
                        "orbital_gap_eV": float(virt_e) - float(occ_e),
                        "weight_pct": float(weight),
                    }
    return transitions


def main():
    parser = argparse.ArgumentParser(
        description="Compare the N lowest excited states across multiple escf runs (TD-DFT and/or BSE/TDA).")
    parser.add_argument("--run", action="append", required=True,
                         metavar="LABEL:PATH",
                         help="run label and path to its escf.out, e.g. hse06-bse-tda:path/escf_bse.out. Repeat.")
    parser.add_argument("--nstates", type=int, default=10,
                         help="number of lowest states to keep per run (default: 10)")
    parser.add_argument("--output", default="low_states_comparison.csv")
    args = parser.parse_args()

    fieldnames = ["run", "state", "energy_eV", "energy_nm", "osc_vel", "osc_len",
                  "occ_orbital", "occ_energy_eV", "virt_orbital", "virt_energy_eV",
                  "orbital_gap_eV", "weight_pct", "binding_energy_eV"]

    all_rows = []
    for entry in args.run:
        label, path = entry.split(":", 1)
        rows = parse_escf_table(path)
        transitions = parse_dominant_transitions(path)
        for row in rows[:args.nstates]:
            merged = {"run": label, **row, **transitions.get(row["state"], {})}
            if "orbital_gap_eV" in merged:
                # Exciton binding energy of the dominant occ->virt pair for
                # this state: orbital_gap_eV - energy_eV (see docstring).
                merged["binding_energy_eV"] = merged["orbital_gap_eV"] - merged["energy_eV"]
            all_rows.append(merged)

    with open(args.output, "w") as f:
        f.write(",".join(fieldnames) + "\n")
        for r in all_rows:
            values = []
            for key in fieldnames:
                value = r.get(key, "")
                if isinstance(value, float):
                    value = f"{value:.5f}"
                values.append(str(value))
            f.write(",".join(values) + "\n")
    print(f"Comparison table saved -> {args.output}")

    print(f"\nS1 (state 1) across runs:")
    for entry in args.run:
        label, path = entry.split(":", 1)
        s1 = parse_escf_table(path)[0]
        s1_transition = parse_dominant_transitions(path).get(1, {})
        gap_info = ""
        if s1_transition:
            gap_info = (f", dominant {s1_transition['occ_orbital']}->{s1_transition['virt_orbital']} "
                        f"({s1_transition['weight_pct']:.1f}%), orbital gap = {s1_transition['orbital_gap_eV']:.4f} eV")
        print(f"  [{label}] {s1['energy_eV']:.4f} eV ({s1['energy_nm']:.1f} nm), "
              f"f(vel) = {s1['osc_vel']:.5f}, f(len) = {s1['osc_len']:.5f}{gap_info}")


if __name__ == "__main__":
    main()
