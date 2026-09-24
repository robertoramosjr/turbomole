#!/usr/bin/env python3
"""
compute_exciton_binding.py

Computes the exciton binding energy, E_b = E_gap(HOMO-LUMO) -
E_exciton(S1), for one or more levels of theory.

For plain TD-DFT levels, E_gap is the Kohn-Sham HOMO-LUMO gap of that
same functional's converged SCF (read from the "HOMO-LUMO gap:" line
of its ridft/dscf output) -- gap and S1 must come from the same
functional, never mixed across runs.

For a GW-BSE/TDA level, E_gap is the G0W0 quasiparticle HOMO-LUMO gap
(read from the second data column of qpenergies.dat, NOT the first
Kohn-Sham column) -- this is the physically correct pairing, since GW
corrects the single-particle gap before BSE adds the exciton
correction on top of it.

Both the Kohn-Sham and G0W0 energies of the requested HOMO/LUMO rows
are printed for manual verification, since the two columns in
qpenergies.dat sit right next to each other and are easy to swap by
mistake. Orbital rows are matched by their label in the first column
(e.g. "81a"), never by line position, so a possible header/offset
mismatch cannot silently pick the wrong row.

Usage:
    python ~/scripts/compute_exciton_binding.py \
        --tddft "TD-DFT/HSE06":artepillin_C_d3_hse06/ridft_tight_d3.out:artepillin_C_d3_hse06/escf_d3.out \
        --tddft "TD-DFT/B3LYP":artepillin_C_d3_b3lyp/ridft_tight_b3lyp.out:artepillin_C_d3_b3lyp/escf_b3lyp.out \
        --tddft "TD-DFT/CAM-B3LYP":artepillin_C_d3_camb3lyp/ridft_tight_camb3lyp.out:artepillin_C_d3_camb3lyp/escf_camb3lyp.out \
        --gwbse "GW-BSE/TDA":artepillin_C_d3_hse06_gwbse/qpenergies.dat:artepillin_C_d3_hse06_gwbse/escf_bse.out:81:82 \
        --output exciton_binding_energy.csv
"""

import argparse
import re

HOMO_LUMO_GAP_RE = re.compile(r"^\s*HOMO-LUMO gap:\s*[-\d.]+\s*H\s*=\s*\+?([-\d.]+)\s*eV")
HOMO_LINE_RE = re.compile(r"^\s*HOMO\s*:\s*[-\d.]+\s*H\s*=\s*\+?([-\d.]+)\s*eV")
LUMO_LINE_RE = re.compile(r"^\s*LUMO\s*:\s*[-\d.]+\s*H\s*=\s*\+?([-\d.]+)\s*eV")

EXCITATION_TABLE_ROW_RE = re.compile(
    r"^\s*\|\s*(\d+)\s+a\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
)


def parse_ks_homo_lumo_gap(scf_out_path):
    """Extract HOMO, LUMO and HOMO-LUMO gap (eV) from a ridft/dscf output."""
    homo_ev = lumo_ev = gap_ev = None
    with open(scf_out_path, "r") as f:
        for line in f:
            if homo_ev is None:
                match = HOMO_LINE_RE.match(line)
                if match:
                    homo_ev = float(match.group(1))
                    continue
            if lumo_ev is None:
                match = LUMO_LINE_RE.match(line)
                if match:
                    lumo_ev = float(match.group(1))
                    continue
            match = HOMO_LUMO_GAP_RE.match(line)
            if match:
                gap_ev = float(match.group(1))
                break
    if gap_ev is None:
        raise ValueError(f"No 'HOMO-LUMO gap:' line found in {scf_out_path}")
    return homo_ev, lumo_ev, gap_ev


def parse_s1_energy(escf_out_path):
    """Extract the state-1 excitation energy (eV) from the SUMMARY table."""
    with open(escf_out_path, "r") as f:
        for line in f:
            match = EXCITATION_TABLE_ROW_RE.match(line)
            if match and int(match.group(1)) == 1:
                return float(match.group(3))
    raise ValueError(f"No state-1 row found in the excitation summary table of {escf_out_path}")


def parse_qp_orbital(qpenergies_path, orbital_index):
    """Return (KS_eV, G0W0_eV) for the row labeled '<orbital_index>a' in qpenergies.dat.

    Matched by the orbital label in the first column, not by line
    position, so an unexpected header/offset cannot pick the wrong row.
    """
    label = f"{orbital_index}a"
    row_re = re.compile(rf"^\s*{orbital_index}a\s+([-\d.]+)\s+([-\d.]+)")
    with open(qpenergies_path, "r") as f:
        for line in f:
            match = row_re.match(line)
            if match:
                return float(match.group(1)), float(match.group(2))
    raise ValueError(f"Orbital '{label}' not found in {qpenergies_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute the exciton binding energy (gap - S1) for one or more levels of theory.")
    parser.add_argument("--tddft", action="append", default=[],
                         metavar="LABEL:SCF_OUT:ESCF_OUT",
                         help="Plain TD-DFT level: label, SCF output with the HOMO-LUMO gap of that "
                              "same functional, and its escf.out with the S1 excitation energy. Repeat.")
    parser.add_argument("--gwbse", action="append", default=[],
                         metavar="LABEL:QPENERGIES:ESCF_OUT:HOMO_INDEX:LUMO_INDEX",
                         help="GW-BSE/TDA level: label, qpenergies.dat, escf.out (BSE) with the S1 "
                              "excitation energy, and the HOMO/LUMO orbital indices to look up in "
                              "qpenergies.dat. Repeat.")
    parser.add_argument("--output", default="exciton_binding_energy.csv")
    args = parser.parse_args()

    if not args.tddft and not args.gwbse:
        parser.error("at least one --tddft or --gwbse entry is required")

    results = []

    for entry in args.tddft:
        label, scf_out, escf_out = entry.split(":", 2)
        homo_ev, lumo_ev, gap_ev = parse_ks_homo_lumo_gap(scf_out)
        s1_ev = parse_s1_energy(escf_out)
        print(f"[{label}] KS HOMO = {homo_ev:.5f} eV, KS LUMO = {lumo_ev:.5f} eV, "
              f"KS gap = {gap_ev:.5f} eV (source: {scf_out})")
        results.append({
            "level_of_theory": label,
            "gap_eV": gap_ev,
            "gap_source": f"KS ({label.split('/')[-1]})",
            "S1_eV": s1_ev,
            "E_binding_eV": gap_ev - s1_ev,
        })

    for entry in args.gwbse:
        label, qpenergies, escf_out, homo_index, lumo_index = entry.split(":", 4)
        homo_ks_ev, homo_g0w0_ev = parse_qp_orbital(qpenergies, homo_index)
        lumo_ks_ev, lumo_g0w0_ev = parse_qp_orbital(qpenergies, lumo_index)
        ks_gap_ev = lumo_ks_ev - homo_ks_ev
        g0w0_gap_ev = lumo_g0w0_ev - homo_g0w0_ev
        print(f"[{label}] orbital {homo_index}a: KS = {homo_ks_ev:.5f} eV, G0W0 = {homo_g0w0_ev:.5f} eV")
        print(f"[{label}] orbital {lumo_index}a: KS = {lumo_ks_ev:.5f} eV, G0W0 = {lumo_g0w0_ev:.5f} eV")
        print(f"[{label}] KS gap = {ks_gap_ev:.5f} eV (NOT used) | G0W0 gap = {g0w0_gap_ev:.5f} eV (used)")
        s1_ev = parse_s1_energy(escf_out)
        results.append({
            "level_of_theory": label,
            "gap_eV": g0w0_gap_ev,
            "gap_source": "G0W0 (quasiparticle)",
            "S1_eV": s1_ev,
            "E_binding_eV": g0w0_gap_ev - s1_ev,
        })

    with open(args.output, "w") as f:
        f.write("level_of_theory,gap_eV,gap_source,S1_eV,E_binding_eV\n")
        for r in results:
            f.write(f"{r['level_of_theory']},{r['gap_eV']:.4f},{r['gap_source']},"
                     f"{r['S1_eV']:.4f},{r['E_binding_eV']:.4f}\n")

    print(f"\nlevel_of_theory,gap_eV,gap_source,S1_eV,E_binding_eV")
    for r in results:
        print(f"{r['level_of_theory']},{r['gap_eV']:.4f},{r['gap_source']},"
              f"{r['S1_eV']:.4f},{r['E_binding_eV']:.4f}")
    print(f"\nSaved -> {args.output}")


if __name__ == "__main__":
    main()
