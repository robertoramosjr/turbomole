#!/usr/bin/env python3
"""Diffuse-function test of the aqueous TD-DFT band: HSE06-D3(BJ)/def2-TZVPD + COSMO ridft + escf.

For the --top conformers of highest w_G of each state (water_gens_conformers.csv, with an existing
water_tddft/<state>/<cNNNN>/ result) it builds <root>/todo_v12/diffuse/<state>/<cNNNN>/ with
  * define template (HSE06, grid m4, senex, D3(BJ), --basis, fresh start orbitals) on the SAME geometry, then
  * the COSMO block (epsilon 78.4, radii, cosmo_correlated) copied from the production water_tddft control, and
  * $scfiterlimit 150.
ridft and escf are run by job_diffuse_tddft.sh (the escf keywords are appended after ridft converged).
A joblist (top-weighted conformer of each state first) is written to todo_v12/diffuse/joblist.txt.

Usage:
  python prepare_diffuse_test.py --root ~/work_turbomole/artepillin_C_studies --top 3 --basis def2-TZVPD
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
from prepare_xtb_gtherm_pilot import coord_to_xyz  # noqa: E402

STATES = {"N": "neutral", "A": "monoanion_a", "B": "monoanion_b", "D": "dianion"}
CHARGE = {"N": 0, "A": -1, "B": -1, "D": -2}


def cosmo_block(text):
    """Lines from `$cosmo` up to (and including) `$cosmo_correlated` of the production control."""
    out, on = [], False
    for ln in text.splitlines():
        if ln.startswith("$cosmo"):
            on = True
        elif ln.startswith("$") and on and not ln.startswith("$cosmo"):
            on = False
        if on:
            out.append(ln)
    return "\n".join(out) + "\n"


def add_jbas(ctl, jbas="def2-TZVP"):
    """define has no RI-J fitting set for def2-TZVPD: attach the def2-TZVP one (as in production) per element."""
    ctl = re.sub(r"^(\s+basis =(\w+) \S+)\s*$", lambda m: f"{m.group(1)}  \\\n   jbas  ={m.group(2)} {jbas}", ctl, flags=re.M)
    return ctl.replace("$end", "$jbas file=auxbasis\n$end", 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--top", type=int, default=3)
    ap.add_argument("--basis", default="def2-TZVPD")
    ap.add_argument("--roots", type=int, default=30, help="number of singlet roots ($soes a N)")
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.root))
    gens = list(csv.DictReader(open(os.path.join(root, "paper_2/Figures/data/water_gens_conformers.csv"))))
    tmpl = os.path.join(HERE, "scripts", "setup_hse06_tzvp_template.py")
    jobs = {st: [] for st in STATES}
    for st, d in STATES.items():
        rows = sorted((r for r in gens if r["state"] == st
                       and os.path.isfile(os.path.join(root, "water_tddft", d, r["job"], "exspectrum"))),
                      key=lambda r: -float(r["w_G"]))[:a.top]
        for r in rows:
            src = os.path.join(root, "water_tddft", d, r["job"])
            dst = os.path.join(root, "todo_v12", "diffuse", st, r["job"])
            if os.path.isfile(os.path.join(dst, "control")):
                jobs[st].append(dst)
                continue
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            xyz = dst + ".xyz"
            open(xyz, "w").write(coord_to_xyz(os.path.join(src, "coord")))
            subprocess.run([sys.executable, tmpl, "--xyz", xyz, "--charge", str(CHARGE[st]), "--basis", a.basis, "--nstates", str(a.roots),
                            "--output-dir", dst], check=True, stdout=subprocess.DEVNULL)
            # the production geometry must be kept bit-for-bit (xyz round trip rewrites coord)
            shutil.copy(os.path.join(src, "coord"), os.path.join(dst, "coord"))
            ctl = open(os.path.join(dst, "control")).read()
            ctl = re.sub(r"^\$scfiterlimit.*\n", "", ctl, flags=re.M)
            blk = cosmo_block(open(os.path.join(src, "control")).read())
            ctl = add_jbas(ctl)
            shutil.copy(os.path.join(src, "auxbasis"), os.path.join(dst, "auxbasis"))
            ctl = blk + ctl.replace("$end", "$scfiterlimit 150\n$end", 1)
            open(os.path.join(dst, "control"), "w").write(ctl)
            print("prepared", st, r["job"], "w_G=%s" % r["w_G"], flush=True)
            jobs[st].append(dst)
    order = [jobs[s][i] for i in range(a.top) for s in STATES if i < len(jobs[s])]
    open(os.path.join(root, "todo_v12", "diffuse", "joblist.txt"), "w").write("\n".join(order) + "\n")
    print(len(order), "directories in joblist")


if __name__ == "__main__":
    main()
