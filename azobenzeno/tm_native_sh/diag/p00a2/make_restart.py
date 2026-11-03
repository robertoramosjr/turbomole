#!/usr/bin/env python3
"""Build a restart directory for the TRAJ2 cycle-256 snapshot (P00a), with a COMPLETE control.

The frog write at cycle 256 truncates `control` (31 lines); the last complete control is the
`updatelog` (manual M1: SH needs $surface_hopping, $currentstate, $nacme, weight derivatives).
`full` keeps all per-cycle blocks of $nac/$excitation; `pruned` keeps the last K of each.

Usage: make_restart.py SNAPSHOT_DIR OUT_DIR [--prune-to K] [--nsteps N]
"""
import argparse
import re
import shutil
from pathlib import Path

STALE = ("GEO_OPT_FAILED", "GEO_OPT_RUNNING", "abend.ridft", "not.converged", "stop", "STOP",
         "cycled_write_complex.tmp", "cycled_write_int.tmp", "cycled_write_real.tmp")
PER_CYCLE_GROUPS = ("$excitation", "$nac ")


def prune_groups(lines: list[str], keep: int) -> list[str]:
    """Keep only the last `keep` `cycle` blocks inside $excitation and $nac."""
    out, i = [], 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(PER_CYCLE_GROUPS):
            j = i + 1
            while j < len(lines) and not lines[j].startswith("$"):
                j += 1
            body = lines[i + 1:j]
            starts = [k for k, b in enumerate(body) if re.match(r"\s*cycle\s*=", b)]
            out.append(line)
            out.extend(body[starts[-keep]:] if len(starts) > keep else body)
            i = j
        else:
            out.append(line)
            i += 1
    return out


def make_restart(snapshot: Path, out: Path, prune_to: int | None = None, nsteps: int = 4) -> None:
    if out.exists():
        raise FileExistsError(out)
    shutil.copytree(snapshot, out, ignore=shutil.ignore_patterns("s0", "NOTAS.md", "job.*", "tmp.output"))
    for name in STALE:
        (out / name).unlink(missing_ok=True)
    lines = (snapshot / "updatelog").read_text().split("\n")
    if prune_to is not None:
        lines = prune_groups(lines, prune_to)
    (out / "control").write_text("\n".join(lines))
    md = (out / "mdmaster").read_text()
    md = re.sub(r"(?m)^\$nsteps\s+\d+", f"$nsteps      {nsteps}", md)
    (out / "mdmaster").write_text(md)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("snapshot", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--prune-to", type=int, default=None, metavar="K")
    ap.add_argument("--nsteps", type=int, default=4)
    a = ap.parse_args()
    make_restart(a.snapshot, a.out, a.prune_to, a.nsteps)


if __name__ == "__main__":
    main()
