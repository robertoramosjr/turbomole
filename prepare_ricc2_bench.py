#!/usr/bin/env python3
"""Excited-state benchmark of the aqueous COSMO-TD-DFT band against RI-ADC(2) and RI-CC2 (post-SCF COSMO, LR).

For the top-weighted aqueous conformer of each state (water_gens_conformers.csv with a water_tddft result) it builds,
on the SAME production geometry (water_tddft/<state>/<cNNNN>/coord), under <root>/todo_v12/adc2/<state>/<cNNNN>/:
  hf/    HF/<basis> + COSMO(eps 78.4, refind 1.3 default) reference for ricc2 (ridft), with $freeze, $ricc2 and
         $excitations (nexc singlets, oscillator strengths through `spectrum` on its own line (it is ignored on the irrep line))
         and `$reaction_field post-SCF ccs-like` (COSMO-ADC(2)/CC2, linear-response fast term, as in the production
         fast-response TD-DFT); the model line ($ricc2 adc(2) or cc2) is written by job_ricc2_bench.sh
  dft/   HSE06-D3(BJ)/<basis> + COSMO TD-DFT (escf, nexc roots) on the same geometry: same-basis reference.
Joblist: todo_v12/adc2/joblist.txt (one hf dir per line; dft is the sibling directory).

Usage:
  python prepare_ricc2_bench.py --root ~/work_turbomole/artepillin_C_studies --basis def2-SVP --nexc 8 [--states N A B D]
"""
import argparse
import csv
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from prepare_diffuse_test import CHARGE, STATES, cosmo_block  # noqa: E402
from prepare_xtb_gtherm_pilot import coord_to_xyz  # noqa: E402

DROP_HF = re.compile(r"^\$(dft|senex|disp3|scfinstab|soes|rpacor|xctype)\b")


def strip_groups(ctl, pattern):
    out, skip = [], False
    for ln in ctl.splitlines():
        if ln.startswith("$"):
            skip = bool(pattern.match(ln))
        if not skip:
            out.append(ln)
    return "\n".join(out) + "\n"


def n_core(coord):
    """frozen-core orbitals = number of heavy atoms (C, N, O: 1s)."""
    return sum(1 for ln in open(coord) if re.search(r"\s[cno]\s*$", ln.rstrip().lower()))


def build(tmpl, xyz, charge, basis, nstates, dst):
    subprocess.run([sys.executable, tmpl, "--xyz", xyz, "--charge", str(charge), "--basis", basis,
                    "--nstates", str(nstates), "--output-dir", dst], check=True, stdout=subprocess.DEVNULL)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--basis", default="def2-SVP")
    ap.add_argument("--nexc", type=int, default=8)
    ap.add_argument("--states", nargs="+", default=list(STATES))
    ap.add_argument("--maxcor", type=int, default=3000, help="MiB per core for $maxcor")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    gens = list(csv.DictReader(open(os.path.join(root, "paper_2/Figures/data/water_gens_conformers.csv"))))
    tmpl = os.path.join(HERE, "scripts", "setup_hse06_tzvp_template.py")
    jobs = []
    for st in a.states:
        d = STATES[st]
        r = max((r for r in gens if r["state"] == st
                 and os.path.isfile(os.path.join(root, "water_tddft", d, r["job"], "exspectrum"))),
                key=lambda r: float(r["w_G"]))
        src = os.path.join(root, "water_tddft", d, r["job"])
        base = os.path.join(root, "todo_v12", "adc2", st, r["job"])
        if os.path.isfile(os.path.join(base, "hf", "control")):
            jobs.append(os.path.join(base, "hf"))
            continue
        os.makedirs(base, exist_ok=True)
        xyz = os.path.join(base, "geom.xyz")
        open(xyz, "w").write(coord_to_xyz(os.path.join(src, "coord")))
        blk = cosmo_block(open(os.path.join(src, "control")).read())
        for kind in ("hf", "dft"):
            dst = os.path.join(base, kind)
            build(tmpl, xyz, CHARGE[st], a.basis, a.nexc if kind == "hf" else 30, dst)
            shutil.copy(os.path.join(src, "coord"), os.path.join(dst, "coord"))   # bit-for-bit production geometry
            ctl = open(os.path.join(dst, "control")).read()
            ctl = re.sub(r"^\$scfiterlimit.*\n", "", ctl, flags=re.M)
            if kind == "hf":
                ctl = strip_groups(ctl, DROP_HF)
                nc = n_core(os.path.join(dst, "coord"))
                extra = (f"$freeze\n implicit core= {nc} virt= 0\n$maxcor {a.maxcor} MiB per_core\n"
                         f"$excitations\n irrep=a  multiplicity=1  nexc={a.nexc}  npre={a.nexc + 1}  nstart={a.nexc + 3}\n"
                         " spectrum  states=all  operators=diplen\n"
                         "$reaction_field post-SCF ccs-like\n")
                ctl = ctl.replace("$end", extra + "$scfiterlimit 150\n$end", 1)
            else:
                ctl = ctl.replace("$end", "$scfiterlimit 150\n$end", 1)
            open(os.path.join(dst, "control"), "w").write(blk + ctl)
        print("prepared", st, r["job"], "w_G=%s" % r["w_G"], "core orbitals", nc, flush=True)
        jobs.append(os.path.join(base, "hf"))
    open(os.path.join(root, "todo_v12", "adc2", "joblist.txt"), "w").write("\n".join(jobs) + "\n")
    print(len(jobs), "benchmark sets in joblist")


if __name__ == "__main__":
    main()
