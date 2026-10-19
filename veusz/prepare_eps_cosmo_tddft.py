"""Prepare COSMO-TDDFT (HSE06/def2-TZVP, eps=78.4, fast-response escf) copies of the top-weighted conformers.

Selection per state: HSE06-labelled conformers by r2scan3c_population until >= --cover of the labelled
population, at most --max-n (--n-fixed overrides).  Each job is copied to
<root>/eps_cosmo_tddft/<N|A>/job_NNNN (originals untouched): from layer_b/cosmo_jobs when a COSMO ridft
exists there, otherwise from the vacuum job with the $cosmo block added.  Writes joblist.txt (absolute paths).
Usage: LC_ALL=C python prepare_eps_cosmo_tddft.py [--n-fixed A=8] [--dry-run]
"""
import argparse, os, re, shutil, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from anion_band_coverage import SDIR, hse_jobs, read_pop  # noqa: E402
from anion_band_states import escf_path  # noqa: E402

KEEP = ("control", "coord", "basis", "auxbasis", "mos")


def select(pop, lab, cover, max_n):
    tot = sum(pop[j] for j in lab); acc, sel = 0.0, []
    for j in sorted(lab, key=lambda j: -pop[j]):
        sel.append(j); acc += pop[j] / tot
        if acc >= cover or len(sel) == max_n:
            break
    return sel


def cosmo_block(coord_path):
    """$cosmo + $cosmo_atoms text (same radii as the production COSMO jobs), elements by coord order."""
    els = []
    for ln in open(coord_path):
        t = ln.split()
        if len(t) == 4 and t[3].isalpha():
            els.append(t[3].lower())
    rad = {"o": 1.72, "c": 2.00, "h": 1.30}
    out = ["$cosmo", "  epsilon=   78.400", "  rsolv= 1.30", "$cosmo_atoms", "# radii in Angstrom units"]
    i = 0
    while i < len(els):
        k = i
        while k + 1 < len(els) and els[k + 1] == els[i]:
            k += 1
        out += [f"{els[i]}  {i + 1}-{k + 1}".ljust(79) + "\\", f"   radius=  {rad[els[i]]:.4f}"]
        i = k + 1
    return "\n".join(out + ["$cosmo_out file=out.ccf", "$cosmo_data file=cosmo_transfer.tmp"]) + "\n"


def prepare(root, state, name, sel, dry):
    d = os.path.join(root, f"artepillin_c_{name}_d3_hse06")
    paths = []
    for j in sel:
        dst = os.path.join(root, "eps_cosmo_tddft", state, j)
        src_c = os.path.join(d, "layer_b", "cosmo_jobs", j)
        src = src_c if os.path.isdir(src_c) else os.path.dirname(escf_path(d, j))
        paths.append(dst)
        print(state, j, "from", "cosmo_jobs" if src == src_c else "vacuum job + $cosmo")
        if dry:
            continue
        os.makedirs(dst, exist_ok=True)
        for f in KEEP:
            shutil.copy(os.path.join(src, f), dst)
        ctl = open(os.path.join(dst, "control")).read()
        if not re.search(r"^\$cosmo\s*$", ctl, re.M):
            ctl = ctl.replace("$end", cosmo_block(os.path.join(dst, "coord")) + "$end")
        if not re.search(r"^\$scfinstab rpas", ctl, re.M) or not re.search(r"^\$soes", ctl, re.M):
            raise SystemExit(f"{dst}: control lacks $scfinstab rpas/$soes")
        open(os.path.join(dst, "control"), "w").write(ctl)
    return paths


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    ap.add_argument("--cover", type=float, default=0.8)
    ap.add_argument("--max-n", type=int, default=10)
    ap.add_argument("--n-fixed", nargs="*", default=[], help="e.g. A=8")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    fixed = dict(x.split("=") for x in a.n_fixed)
    allp = []
    for s, name in SDIR.items():
        d = os.path.join(a.root, f"artepillin_c_{name}_d3_hse06")
        pop = read_pop(d); lab = hse_jobs(d) & set(pop)
        sel = select(pop, lab, a.cover, int(fixed[s]) if s in fixed else a.max_n)
        if s in fixed:
            sel = sorted(lab, key=lambda j: -pop[j])[: int(fixed[s])]
        allp += prepare(a.root, s, name, sel, a.dry_run)
    if not a.dry_run:
        with open(os.path.join(a.root, "eps_cosmo_tddft", "joblist.txt"), "w") as f:
            f.write("\n".join(allp) + "\n")
    print(len(allp), "jobs")


if __name__ == "__main__":
    main()
