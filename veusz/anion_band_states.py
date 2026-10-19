"""H2: character of the lowest/bright HSE06 roots of the top-weighted conformers of N and A.

Writes Figures/data/anion_band_states.csv: for each conformer, the 3 lowest roots plus roots with
f>0.05 in [3.5,5.0] eV, with the dominant orbital pairs (occupied/virtual orbital energies, %).
A positive virtual-orbital energy flags a vacuum-unbound (diffuse) acceptor orbital.
Run with: LC_ALL=C ~/.conda/envs/veusz_env/bin/python anion_band_states.py
"""
import argparse, csv, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from anion_band_coverage import SDIR, hse_jobs, read_pop  # noqa: E402

CONTRIB = re.compile(r"^\s*(\d+)\s+a\s+([-\d.]+)\s+(\d+)\s+a\s+([-\d.]+)\s+([\d.]+)\s*$")


def parse_roots(path):
    roots, cur, in_dom = [], None, False
    for ln in open(path):
        if "singlet a excitation" in ln:
            cur = dict(root=int(ln.split()[0]), f=None, e=None, c=[]); roots.append(cur); in_dom = False
        elif cur is None:
            continue
        elif ln.startswith(" Excitation energy / eV"):
            cur["e"] = float(ln.split(":")[1])
        elif "length representation:" in ln and cur["f"] is None:
            cur["f"] = float(ln.split(":")[1])
        elif "Dominant contributions" in ln:
            in_dom = True
        elif in_dom:
            m = CONTRIB.match(ln)
            if m:
                cur["c"].append((int(m.group(1)), float(m.group(2)), int(m.group(3)), float(m.group(4)), float(m.group(5))))
            elif ln.strip() == "" and cur["c"]:
                in_dom = False
    return roots


def escf_path(d, j):
    for sb in ("labeled_jobs", "confirmation_batch_jobs", "random_control_batch_jobs"):
        p = os.path.join(d, "layer_b", sb, j, "escf.out")
        if os.path.isfile(p):
            return p


def select(roots):
    return [r for r in roots if r["root"] <= 3 or (r["f"] > 0.05 and 3.5 <= r["e"] <= 5.0)]


def run(root, ntop):
    rows = []
    for s, name in SDIR.items():
        d = os.path.join(root, f"artepillin_c_{name}_d3_hse06")
        pop = read_pop(d)
        for j in sorted(hse_jobs(d) & set(pop), key=lambda j: -pop[j])[:ntop]:
            for r in select(parse_roots(escf_path(d, j))):
                for k, (o, eo, v, ev, pc) in enumerate(r["c"][:3]):
                    rows.append(dict(state=s, job=j, weight=round(pop[j], 4), root=r["root"], e_eV=round(r["e"], 3),
                                     f=round(r["f"], 4), occ=o, occ_eV=eo, virt=v, virt_eV=ev, pct=pc, rank=k + 1))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--root", default=os.path.abspath(os.path.join(here, "..", "..")))
    ap.add_argument("--ntop", type=int, default=5)
    ap.add_argument("--out", default=os.path.join(here, "..", "Figures", "data", "anion_band_states.csv"))
    a = ap.parse_args()
    rows = run(a.root, a.ntop)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(len(rows), "rows ->", a.out)


if __name__ == "__main__":
    main()
