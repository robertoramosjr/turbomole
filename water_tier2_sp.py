#!/usr/bin/env python3
"""
water_tier2_sp.py

Tier 2 of the solvent-aware ensemble (decisive filter): r2SCAN-3c/def2-mTZVP + COSMO(water) SINGLE POINTS on every
conformer of a GFN2/ALPB CREST ensemble, replacing the GFN2 ranking (which mis-ranks charged conformers) before any
population cut. Two sub-commands:

  template : build a COSMO-ready r2SCAN-3c template (control/basis/auxbasis/mos + $cosmo) from an existing converged
             r2SCAN-3c conformer directory of the same protonation state (cosmoprep via setup_cosmo_job.py).
  run      : process one chunk of the ensemble: for each conformer k (1-based, k % nchunks == chunk) write coord from
             crest_conformers.xyz, run ridft, store the total energy (incl. COSMO) in <outdir>/res_chunkNN.csv.

Usage:
  python water_tier2_sp.py template --source layer_a/conformer_jobs/job_0044 --out tier2/neutral/template
  python water_tier2_sp.py run --xyz water_ensemble/neutral/crest_conformers.xyz --template tier2/neutral/template \
         --outdir tier2/neutral --chunk 3 --nchunks 40
"""
import argparse
import csv
import os
import re
import shutil
import subprocess

BOHR = 1.0 / 0.529177210903


def cmd_template(a):
    sc = os.path.join(os.path.dirname(os.path.abspath(__file__)), "setup_cosmo_job.py")
    for attempt in range(3):
        shutil.rmtree(a.out, ignore_errors=True)
        subprocess.run(["python3", sc, "--source-dir", a.source, "--output-dir", a.out, "--epsilon", "78.4"], check=False)
        c = os.path.join(a.out, "control")
        if os.path.isfile(c) and re.search(r"^\$cosmo$", open(c).read(), re.M) and "amat file=*" not in open(c).read():
            print("template OK", a.out)
            return
    raise SystemExit("template failed: " + a.out)


def read_xyz_frames(path):
    lines = open(path).read().split("\n")
    nat = int(lines[0])
    frames = []
    for i in range(0, len(lines) - nat, nat + 2):
        if not lines[i].strip():
            break
        frames.append([l.split() for l in lines[i + 2:i + 2 + nat]])
    return frames


def cmd_run(a):
    frames = read_xyz_frames(a.xyz)
    os.makedirs(a.outdir, exist_ok=True)
    res = os.path.join(a.outdir, f"res_chunk{a.chunk:02d}.csv")
    done = set()
    if os.path.exists(res):
        done = {int(r["crest_index"]) for r in csv.DictReader(open(res))}
    out = open(res, "a", newline="")
    w = csv.writer(out)
    if not done:
        w.writerow(["crest_index", "E_total_Eh", "converged", "iterations"])
    for k in range(1, len(frames) + 1):
        if (k - 1) % a.nchunks != a.chunk - 1 or k in done:
            continue
        d = os.path.join(a.outdir, f"c{k:04d}")
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        for f in ("control", "basis", "auxbasis", "mos"):
            shutil.copy(os.path.join(a.template, f), d)
        ctrl = open(os.path.join(d, "control")).read()
        ctrl = re.sub(r"^\$scfiterlimit.*\n", "", ctrl, flags=re.M)
        ctrl = ctrl.replace("$end", "$scfiterlimit 300\n$end", 1) if "$scfiterlimit" not in ctrl else ctrl
        # drop stale optimisation bookkeeping that is irrelevant for a single point
        ctrl = re.sub(r"^\$(last step|actual step).*\n", "", ctrl, flags=re.M)
        open(os.path.join(d, "control"), "w").write(ctrl)
        with open(os.path.join(d, "coord"), "w") as f:
            f.write("$coord\n")
            for at in frames[k - 1]:
                f.write(f"{float(at[1]) * BOHR:20.12f}{float(at[2]) * BOHR:20.12f}{float(at[3]) * BOHR:20.12f}     {at[0].lower()}\n")
            f.write("$end\n")
        r = subprocess.run("ridft", cwd=d, capture_output=True, text=True).stdout
        e = re.findall(r"total energy\s*=\s*(-?\d+\.\d+)", r)
        conv = "convergence criteria satisfied" in r
        nit = len(re.findall(r"^\s*ITERATION", r, re.M))
        w.writerow([k, e[-1] if e else "", int(conv), nit])
        out.flush()
        shutil.rmtree(d, ignore_errors=True)
    out.close()


def main():
    ap = argparse.ArgumentParser(description="Tier-2 r2SCAN-3c/COSMO single points on a CREST ensemble.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("template")
    t.add_argument("--source", required=True)
    t.add_argument("--out", required=True)
    r = sub.add_parser("run")
    r.add_argument("--xyz", required=True)
    r.add_argument("--template", required=True)
    r.add_argument("--outdir", required=True)
    r.add_argument("--chunk", type=int, required=True)
    r.add_argument("--nchunks", type=int, required=True)
    a = ap.parse_args()
    (cmd_template if a.cmd == "template" else cmd_run)(a)


if __name__ == "__main__":
    main()
