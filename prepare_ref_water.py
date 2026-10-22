#!/usr/bin/env python3
"""Re-run the pKa reference compounds with the aqueous-ensemble protocol of Paper 2.

Same chain as the artepillin C conformers (collect_water_gens.py): r2SCAN-3c/def2-mTZVP + COSMO(water)
geometry re-optimisation, then HSE06-D3(BJ)/def2-TZVP + COSMO single point (W4, energy with the outlying-charge
correction) and GFN2-xTB `--hess` G_thermal at that geometry.  Starting geometries are the fully optimised HSE06
structures of pka_reference/<ref>/<acid|base>/coord.

Sub-commands (all light; run from the login node, the heavy steps are sbatch jobs):
  reopt : build <root>/ref_water/<ref>/<acid|base>/reopt (define r2SCAN-3c + cosmoprep) and write reopt_joblist.txt
  w4    : for converged reopt dirs build <ref>/<kind>/w4 (HSE06 COSMO control from pka_reference/<ref>/<kind>_cosmo)
          and xtb/ (mol.xyz, charge.txt); write w4_joblist.txt and xtb_joblist.txt

Usage:
  python prepare_ref_water.py reopt --root ~/work_turbomole/artepillin_C_studies
  python prepare_ref_water.py w4    --root ~/work_turbomole/artepillin_C_studies
"""
import argparse
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from prepare_xtb_gtherm_pilot import coord_to_xyz  # noqa: E402
from prepare_w4_sp import clean_control  # noqa: E402

# reference -> (directory in pka_reference, charge of acid, charge of base)
REFS = {
    "acetic_acid": ("acetic_acid", 0, -1),
    "phenol": ("phenol", 0, -1),
    "26_dimethylphenol": ("26_dimethylphenol", 0, -1),
    "4_hydroxybenzoate": ("4_hydroxybenzoate", -1, -2),
    "4_coumarate": ("4_coumarate", -1, -2),
}


def species(root):
    for ref, (src, qa, qb) in REFS.items():
        for kind, q in (("acid", qa), ("base", qb)):
            yield ref, kind, q, os.path.join(root, "pka_reference", src, kind), os.path.join(root, "ref_water", ref, kind)


def cmd_reopt(a):
    root = os.path.abspath(os.path.expanduser(a.root))
    jobs = []
    for ref, kind, q, src, dst in species(root):
        if os.path.isfile(os.path.join(dst, "reopt", "control")):
            jobs.append(os.path.join(dst, "reopt"))
            continue
        os.makedirs(dst, exist_ok=True)
        xyz = os.path.join(dst, "start.xyz")
        open(xyz, "w").write(coord_to_xyz(os.path.join(src, "coord")))
        tmpl = os.path.join(dst, "template")
        shutil.rmtree(tmpl, ignore_errors=True)
        subprocess.run([sys.executable, os.path.join(HERE, "scripts", "setup_r2scan3c_template.py"), "--xyz", xyz,
                        "--charge", str(q), "--output-dir", tmpl], check=True, stdout=subprocess.DEVNULL)
        for attempt in range(3):
            shutil.rmtree(os.path.join(dst, "reopt"), ignore_errors=True)
            subprocess.run([sys.executable, os.path.join(HERE, "setup_cosmo_job.py"), "--source-dir", tmpl,
                            "--output-dir", os.path.join(dst, "reopt"), "--epsilon", "78.4"], stdout=subprocess.DEVNULL)
            c = os.path.join(dst, "reopt", "control")
            if os.path.isfile(c) and "$cosmo\n" in open(c).read() and "amat file=*" not in open(c).read():
                break
        else:
            raise SystemExit(f"cosmoprep failed for {ref}/{kind}")
        shutil.copy(os.path.join(tmpl, "coord"), os.path.join(dst, "reopt", "coord"))
        print("reopt ready", ref, kind, "charge", q)
        jobs.append(os.path.join(dst, "reopt"))
    open(os.path.join(root, "ref_water", "reopt_joblist.txt"), "w").write("\n".join(jobs) + "\n")


def cmd_w4(a):
    root = os.path.abspath(os.path.expanduser(a.root))
    w4, xt = [], []
    for ref, kind, q, src, dst in species(root):
        ro = os.path.join(dst, "reopt")
        if not os.path.isfile(os.path.join(ro, "GEO_OPT_CONVERGED")):
            print("skip (not converged)", ref, kind)
            continue
        tm = os.path.join(root, "pka_reference", REFS[ref][0], kind + "_cosmo")
        d = os.path.join(dst, "w4")
        os.makedirs(d, exist_ok=True)
        for f in ("basis", "auxbasis", "mos"):
            shutil.copy(os.path.join(tm, f), d)
        open(os.path.join(d, "control"), "w").write(clean_control(open(os.path.join(tm, "control")).read(), 100))
        shutil.copy(os.path.join(ro, "coord"), d)
        x = os.path.join(dst, "xtb")
        os.makedirs(x, exist_ok=True)
        open(os.path.join(x, "mol.xyz"), "w").write(coord_to_xyz(os.path.join(ro, "coord")))
        open(os.path.join(x, "charge.txt"), "w").write(f"{q}\n")
        w4.append(d)
        xt.append(x)
    open(os.path.join(root, "ref_water", "w4_joblist.txt"), "w").write("\n".join(w4) + "\n")
    open(os.path.join(root, "ref_water", "xtb_joblist.txt"), "w").write("\n".join(xt) + "\n")
    print(len(w4), "species prepared for W4/xtb")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=("reopt", "w4"))
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    {"reopt": cmd_reopt, "w4": cmd_w4}[a.step](a)


if __name__ == "__main__":
    main()
