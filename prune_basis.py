#!/usr/bin/env python3
"""prune_basis.py -- expand a TURBOMOLE library basis (resolving '->' references) and drop diffuse
uncontracted primitives, writing a local 'basis' file for periodic (riper) runs.

Diffuse functions cause near-linear dependencies and BSSE in dense solids; the usual practice is to
remove uncontracted primitives with exponent below a cutoff (default 0.1 bohr^-2). Contracted shells
are kept whole. Optionally copies ECP blocks into a $ecp group of the same file.

Usage:
  prune_basis.py --elements ag se o --basis def2-TZVP --cutoff 0.1 --suffix p01 \
                 --ecp ag:def2-ecp --output basis
  Then in control: '$basis file=basis', '$ecp file=basis' and in $atoms 'basis = ag def2-TZVP-p01'.

Options:
  --basen DIR    basis library (default: $TURBODIR/basen)
  --report       only print kept/removed shells, do not write
"""
import argparse
import os
import re
import sys
from pathlib import Path


def read_library(path):
    """Return {name: [lines]} for every block in a TURBOMOLE basen file (aliases share the block)."""
    blocks, lines = {}, Path(path).read_text(errors="replace").splitlines()
    i, n = 0, len(lines)
    name_re = re.compile(r"^[a-z]{1,2}\s+\S")
    while i < n:
        if name_re.match(lines[i]):
            names = []
            while i < n and not lines[i].startswith("*"):
                if name_re.match(lines[i]):
                    names.append(" ".join(lines[i].split()))
                i += 1
            i += 1  # opening '*'
            body = []
            while i < n and not lines[i].strip().startswith("*"):
                body.append(lines[i])
                i += 1
            for nm in names:
                blocks.setdefault(nm.lower(), body)
        i += 1
    return blocks


def expand(name, lib, depth=0):
    """Resolve '->' references recursively; return list of (l, [(exp, coef)])."""
    if depth > 10:
        sys.exit(f"reference loop at {name}")
    body = lib.get(name.lower())
    if body is None:
        sys.exit(f"basis '{name}' not found in library")
    shells, k = [], 0
    while k < len(body):
        s = body[k].split()
        if not s or s[0].startswith("#"):
            k += 1
            continue
        if s[0] == "->":
            shells += expand(" ".join(s[1:3]), lib, depth + 1)
            k += 1
            continue
        nprim, l = int(s[0]), s[1].lower()
        prims = [tuple(float(x.replace("D", "E").replace("d", "e")) for x in body[k + 1 + j].split()[:2])
                 for j in range(nprim)]
        shells.append((l, prims))
        k += 1 + nprim
    return shells


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--elements", nargs="+", required=True)
    p.add_argument("--basis", required=True, help="library basis name, e.g. def2-TZVP")
    p.add_argument("--cutoff", type=float, default=0.1, help="drop uncontracted primitives below this exponent")
    p.add_argument("--suffix", default="p01", help="new name = <basis>-<suffix>")
    p.add_argument("--ecp", nargs="*", default=[], help="el:ecpname to copy into $ecp, e.g. ag:def2-ecp")
    p.add_argument("--basen", default=os.path.join(os.environ.get("TURBODIR", ""), "basen"))
    p.add_argument("--output", default="basis")
    p.add_argument("--report", action="store_true")
    p.add_argument("--force", action="store_true", help="overwrite output")
    a = p.parse_args()

    out = ["$basis"]
    for el in a.elements:
        lib = read_library(Path(a.basen) / el)
        shells = expand(f"{el} {a.basis}", lib)
        keep = [sh for sh in shells if not (len(sh[1]) == 1 and sh[1][0][0] < a.cutoff)]
        drop = [sh for sh in shells if sh not in keep]
        cnt = lambda ss: "".join(f"{sum(1 for x in ss if x[0] == l)}{l}" for l in "spdfg" if any(x[0] == l for x in ss))
        print(f"{el}: {cnt(shells)} -> {cnt(keep)}; removed " +
              (", ".join(f"{l}({pr[0][0]:.4g})" for l, pr in drop) or "none"))
        out += ["*", f"{el} {a.basis}-{a.suffix}", f"# {a.basis} without uncontracted exponents < {a.cutoff}", "*"]
        for l, prims in keep:
            out.append(f"   {len(prims)}  {l}")
            out += [f"  {e:20.10f}  {c:20.10f}" for e, c in prims]
    out.append("*")
    if a.ecp:
        out.append("$ecp")
        for spec in a.ecp:
            el, nm = spec.split(":")
            body = read_library(Path(a.basen) / el).get(f"{el} {nm}".lower())
            if body is None:
                sys.exit(f"ecp {spec} not found")
            out += ["*", f"{el} {nm}", "*"] + [l for l in body if not l.startswith("#")]
        out.append("*")
    out.append("$end")
    if a.report:
        return
    dest = Path(a.output)
    if dest.exists() and not a.force:
        sys.exit(f"{dest} exists; use --force")
    dest.write_text("\n".join(out) + "\n")
    print("written", dest)


if __name__ == "__main__":
    main()
