"""Run: python3 test_make_restart.py  (pytest is not installed on access; see PROBLEMAS P01)."""
import re
import sys
import tempfile
from pathlib import Path

from make_restart import make_restart

SNAP = Path(__file__).resolve().parent.parent / "trans_step256"


def _build(prune_to):
    d = Path(tempfile.mkdtemp(dir=Path(__file__).resolve().parent)) / "r"
    make_restart(SNAP, d, prune_to)
    return d, (d / "control").read_text()


def _n_cycle_blocks(text, group):
    m = re.search(rf"(?ms)^\{group}\b.*?(?=^\$[a-z])", text)
    return len(re.findall(r"(?m)^\s*cycle\s*=", m.group(0)))


def test_control_has_sh_groups():
    for prune in (None, 2):
        _, c = _build(prune)
        for key in (r"^\$surface_hopping", r"^\$currentstate", r"^\$nacme", r"weight derivatives", r"^\$atoms"):
            assert len(re.findall(key, c, re.M)) == 1, (prune, key)


def test_block_counts():
    _, full = _build(None)
    _, pruned = _build(2)
    assert _n_cycle_blocks(full, "$nac") == 256 and _n_cycle_blocks(full, "$excitation") == 256
    assert _n_cycle_blocks(pruned, "$nac") == 2 and _n_cycle_blocks(pruned, "$excitation") == 2


def test_stale_and_nsteps():
    d, _ = _build(None)
    assert not (d / "GEO_OPT_FAILED").exists() and not (d / "abend.ridft").exists()
    assert re.search(r"(?m)^\$nsteps\s+4\b", (d / "mdmaster").read_text())
    assert (d / "control").read_text().rstrip().endswith("$end")


if __name__ == "__main__":
    fails = 0
    for n, f in list(globals().items()):
        if n.startswith("test_"):
            try:
                f(); print("ok  ", n)
            except AssertionError as e:
                fails += 1; print("FAIL", n, e)
    sys.exit(fails)
