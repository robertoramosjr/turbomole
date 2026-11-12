#!/usr/bin/env python3
"""gen_riper.py -- build a TURBOMOLE riper run dir (hand-written control, no define) for Ag2SeO3.
Usage: gen_riper.py --geom geom/exp --name pbesol_sp --func pbesol|pbe|b3-lyp --kmesh 5 3 4 --mode sp|opt
                    [--sigma-ev 0.01] [--bands FILE] [--maxcor-mb 500]
Basis pob-TZVP-rev2 (Ag + pob-ecp-rev2 ECP, jbas universal-ecp-28; Se,O all-electron, jbas universal).
sigma is given in eV and written in hartree. PBEsol = LibXC 116 (X) + 133 (C). --bands: file with kptlines/recipr
lines that replaces 'nkpoints' (do the band run as a separate SP after the SCF mesh run)."""
import argparse, re, shutil
from pathlib import Path
EV = 1 / 27.211386245988
p = argparse.ArgumentParser()
p.add_argument("--geom", required=True); p.add_argument("--name", required=True)
p.add_argument("--func", required=True, choices=["pbesol", "pbe", "b3-lyp"])
p.add_argument("--kmesh", type=int, nargs=3, required=True); p.add_argument("--mode", choices=["sp", "opt"], required=True)
p.add_argument("--sigma-ev", type=float, default=0.01); p.add_argument("--bands")
p.add_argument("--maxcor-mb", type=int, default=500)
a = p.parse_args()
d = Path("runs") / a.name; d.mkdir(parents=True, exist_ok=True)
src = (Path(a.geom) / "coord").read_text()
cell = " ".join(re.search(r"\$cell angs\s+(.+)", src).group(1).split())
atoms = [l for l in src.splitlines() if len(l.split()) == 4]
syms = [l.split()[3] for l in atoms]
assert syms == ["ag"] * 8 + ["se"] * 4 + ["o"] * 12, "coord must be ordered ag(8) se(4) o(12)"
(d / "coord").write_text("$coord\n" + "\n".join(atoms) + "\n$end\n")
dft = {"pbesol": " functional libxc 116\n functional libxc add 1 133", "pbe": " functional pbe",
       "b3-lyp": " functional b3-lyp"}[a.func]
if a.mode == "opt": dft += "\n weight derivatives"
kp = Path(a.bands).read_text().rstrip() if a.bands else "nkpoints %d %d %d" % tuple(a.kmesh)
lines = [f"$title Ag2SeO3 riper {a.func} {a.mode}", "$operating system unix", "$symmetry c1", "$coord file=coord",
         "$periodic 3", f"$cell angs\n {cell}",
         "$atoms\nag 1-8 \\\n  basis = ag pob-TZVP-rev2 \\\n  ecp   = ag pob-ecp-rev2 \\\n  jbas  = ag universal-ecp-28\n"
         "se 9-12 \\\n  basis = se pob-TZVP-rev2 \\\n  jbas  = se universal\n"
         "o  13-24 \\\n  basis = o pob-TZVP-rev2 \\\n  jbas  = o universal",
         f"$dft\n{dft}", f"$kpoints\n {kp}",
         "$riper\n sigma %.8e" % (a.sigma_ev * EV) + ("\n lenonly on" if a.mode == "sp" else ""),
         "$scfiterlimit 300", "$scfconv 7", f"$maxcor {a.maxcor_mb} MiB per_core"]
if a.mode == "opt": lines.append("$optcell")
(d / "control").write_text("\n".join(lines) + "\n$end\n")
print("written", d / "control")
