#!/usr/bin/env python3
"""
run_freeh_thermal_correction.py

Runs Turbomole's `freeh` on an already-computed Hessian (from B1's
`aoforce` -- no new DFT) to get the thermal correction to the Gibbs
free energy (ZPE + translational/rotational/vibrational partition
function contribution) at T=298.15 K, p=0.1 MPa -- needed for the pKa
thermodynamic cycle (protonation_and_conformer_investigation.md
Sec. 9.32-9.33): G(gas) = E_elec + G_thermal_correction.

`freeh` is interactive (like `define`) with THREE prompts per run,
confirmed by direct exploration in a scratch copy (2026-09-14) --
mapped in this exact order:
  1. Symmetry number confirmation ("Hit RETURN to accept") -- accept
     the default (this project's conformers are all C1, sigma=1).
  2. Frequency scaling factor ("enter new value for corr") -- accept
     the default (BP/SVP scaling, ~0.9914) -- NOT changed here since
     this project has never applied scaling elsewhere either; flagged
     as a known simplification, not a considered choice.
  3. Temperature/pressure range -- accept the default (298.15 K,
     0.1 MPa) via a blank line, THEN quit with "q" on the SAME prompt
     appearing a second time (freeh loops asking for another T/P
     range after printing results).

Works on a COPY of the job directory (default: a temp dir), never
in-place on a production labeled_jobs directory -- confirmed
empirically that running freeh in-place mutates `control` (adds a
"$last step freeh" bookkeeping line); harmless in itself, but there is
no reason to touch production directories for a read-only quantity
extraction.

Usage:
    python ~/work_turbomole/scripts/run_freeh_thermal_correction.py \
        --job-list job_dirs.txt \
        --labeled-jobs layer_b/labeled_jobs \
        --output layer_b/freeh_thermal_correction.csv
"""

import argparse
import csv
import os
import re
import shutil
import tempfile

import pexpect

ZPE_RE = re.compile(r"zpe=\s*([-\d.]+)\s*kJ/mol")
RESULT_ROW_RE = re.compile(
    r"^\s*([\d.]+)\s+([\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$"
)

FREEH_FILES = ("control", "coord", "vibspectrum", "hessian", "vib_normal_modes")


def run_freeh(work_dir, timeout=30):
    child = pexpect.spawn("freeh", cwd=work_dir, timeout=timeout, encoding="utf-8")
    patterns = [
        "enter new value for corr",
        "Hit RETURN to accept",
        "or enter  q  or  \\*  to quit",
        pexpect.EOF, pexpect.TIMEOUT,
    ]
    transcript = ""
    seen_prompt_twice = False
    while True:
        idx = child.expect(patterns)
        transcript += child.before
        if idx <= 1:
            child.sendline("")
        elif idx == 2:
            if seen_prompt_twice:
                child.sendline("q")
            else:
                child.sendline("")
                seen_prompt_twice = True
        else:
            break
    transcript += child.before if child.before else ""
    return transcript


def parse_thermal_correction(transcript):
    zpe_match = ZPE_RE.search(transcript)
    zpe_kj = float(zpe_match.group(1)) if zpe_match else None

    # the RRHO results table: "T  p  ln(qtrans) ln(qrot) ln(qvib) chem.pot. energy entropy"
    chem_pot_kj = entropy_kj = None
    for line in transcript.splitlines():
        m = RESULT_ROW_RE.match(line)
        if m and abs(float(m.group(1)) - 298.15) < 0.5:
            chem_pot_kj = float(m.group(6))
            entropy_kj = float(m.group(8))
            break
    return zpe_kj, chem_pot_kj, entropy_kj


def main():
    parser = argparse.ArgumentParser(
        description="Extract the thermal Gibbs free energy correction "
                    "(freeh, on existing Hessians -- no new DFT) for a set "
                    "of already-labeled conformers.")
    parser.add_argument("--job-list", required=True)
    parser.add_argument("--labeled-jobs", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.job_list) as f:
        job_dirs = [line.strip() for line in f if line.strip()]

    rows = []
    for jd in job_dirs:
        src = os.path.join(args.labeled_jobs, jd)
        if not all(os.path.isfile(os.path.join(src, fn)) for fn in FREEH_FILES):
            print(f"WARNING: {jd} missing one of {FREEH_FILES}, skipping")
            continue
        with tempfile.TemporaryDirectory() as tmp:
            for fn in FREEH_FILES:
                shutil.copy(os.path.join(src, fn), tmp)
            # strip any stale "$last step" marker so freeh doesn't see a
            # different prompt sequence depending on prior state
            control_path = os.path.join(tmp, "control")
            with open(control_path) as f:
                lines = [l for l in f if not l.strip().startswith("$last step")]
            with open(control_path, "w") as f:
                f.writelines(lines)

            transcript = run_freeh(tmp)
        zpe_kj, chem_pot_kj, entropy_kj = parse_thermal_correction(transcript)
        if zpe_kj is None or chem_pot_kj is None:
            print(f"WARNING: could not parse freeh output for {jd}")
            continue
        # chem.pot. IS the thermal Gibbs correction (ZPE - RT*ln(qtrans*qrot*qvib)),
        # in kJ/mol; convert to Hartree for direct addition to electronic energy
        chem_pot_hartree = chem_pot_kj / 2625.499639
        rows.append((jd, zpe_kj, chem_pot_kj, chem_pot_hartree, entropy_kj))
        print(f"{jd}: ZPE={zpe_kj:.2f} kJ/mol, G_thermal_correction="
              f"{chem_pot_kj:.2f} kJ/mol ({chem_pot_hartree:.6f} Ha), "
              f"S={entropy_kj:.5f} kJ/mol/K")

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["job_dir", "zpe_kJ_mol", "g_thermal_correction_kJ_mol",
                          "g_thermal_correction_hartree", "entropy_kJ_mol_K"])
        for row in rows:
            writer.writerow([row[0], f"{row[1]:.4f}", f"{row[2]:.4f}",
                              f"{row[3]:.8f}", f"{row[4]:.6f}"])

    print(f"\nExtracted {len(rows)} / {len(job_dirs)} conformers.")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
