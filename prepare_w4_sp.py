#!/usr/bin/env python3
"""
prepare_w4_sp.py

Build HSE06/def2-TZVP + COSMO(water) single-point dirs (W4) on water-reoptimised geometries (reopt/<state>/cNNNN/coord).
Control/basis/auxbasis/mos come from a converged HSE06 COSMO job of the same protonation state (mos = SCF guess only);
TD-DFT keywords are stripped and the stale bookkeeping removed.  Writes <root>/w4_sp/<tag>/cNNNN and a joblist.

Usage:
  python prepare_w4_sp.py --state neutral --template eps_cosmo_tddft/N/job_0044 --indices 1 5 --root ~/work_turbomole/artepillin_C_studies --joblist w4_sp/pilot_joblist.txt
"""
import argparse
import os
import re
import shutil

DROP = re.compile(r"^\$(scfinstab|soes|rpacor|last|actual|subenergy|charge from|dipole from|energy |grad |dipgrad|hessian|vibrational|nvibro|orbital_max_rnorm|forceapprox)")


def clean_control(text, scfiterlimit):
    out, skip = [], False
    for ln in text.split("\n"):
        if ln.startswith("$"):
            skip = bool(DROP.match(ln)) or ln.startswith(("$scfiterlimit", "$end"))
        if not skip and ln.strip():
            out.append(ln)
    # keep the nuclear-energy/subenergy lines out; add the iteration cap and terminator
    out = [l for l in out if not re.match(r"^\s*-?\d+\.\d+\s+-?\d+\.\d+\s+", l)]
    out.append(f"$scfiterlimit {scfiterlimit}")
    out.append("$end")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--state", required=True, help="reopt subdir, e.g. neutral|monoanion_a|monoanion_b|dianion")
    ap.add_argument("--template", required=True, help="converged HSE06 COSMO job dir (relative to --root)")
    ap.add_argument("--indices", type=int, nargs="+", required=True, help="conformer numbers (reopt/<state>/cNNNN)")
    ap.add_argument("--root", required=True)
    ap.add_argument("--joblist", required=True, help="joblist file (relative to --root), appended to")
    ap.add_argument("--scfiterlimit", type=int, default=100)
    a = ap.parse_args()
    root = os.path.expanduser(a.root)
    tpl = os.path.join(root, a.template)
    jl = []
    for i in a.indices:
        src = os.path.join(root, "reopt", a.state, f"c{i:04d}")
        if not os.path.isfile(os.path.join(src, "GEO_OPT_CONVERGED")):
            print(f"skip c{i:04d}: no GEO_OPT_CONVERGED in {src}")
            continue
        dst = os.path.join(root, "w4_sp", a.state, f"c{i:04d}")
        os.makedirs(dst, exist_ok=True)
        for f in ("basis", "auxbasis", "mos"):
            shutil.copy(os.path.join(tpl, f), dst)
        shutil.copy(os.path.join(src, "coord"), dst)
        with open(os.path.join(tpl, "control")) as fh:
            ctl = clean_control(fh.read(), a.scfiterlimit)
        with open(os.path.join(dst, "control"), "w") as fh:
            fh.write(ctl)
        jl.append(dst)
    with open(os.path.join(root, a.joblist), "a") as fh:
        fh.write("".join(d + "\n" for d in jl))
    print(f"prepared {len(jl)} dirs for {a.state}")


if __name__ == "__main__":
    main()
