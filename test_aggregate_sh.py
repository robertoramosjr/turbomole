import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aggregate_sh import wilson, aggregate


def test_wilson_pilot():            # piloto cis: 3/10 -> 11-60 % (registrado no cofre)
    lo, hi = wilson(3, 10)
    assert round(lo, 2) == 0.11 and round(hi, 2) == 0.60


def test_wilson_empty():
    assert wilson(0, 0) == (0.0, 0.0)


def test_denominator_only_hop_ok():  # Review Focus 5
    rows = [{"status": "hop_ok", "product": "trans"}, {"status": "hop_ok", "product": "cis"},
            {"status": "no_hop", "product": None}, {"status": "died_S1", "product": None}]
    a = aggregate(rows, frm="cis")
    assert a["n_hop_ok"] == 2 and a["frac_switch"] == 0.5
    assert (a["n_total"], a["n_no_hop"], a["n_died_S1"], a["n_died_S0"]) == (4, 1, 1, 0)
    assert (a["n_cis"], a["n_trans"]) == (1, 1)


def test_direction_trans():
    rows = [{"status": "hop_ok", "product": "cis"}] * 3 + [{"status": "hop_ok", "product": "trans"}]
    a = aggregate(rows, frm="trans")
    assert a["frac_switch"] == 0.75


def test_no_hop_ok_gives_none_fraction():
    a = aggregate([{"status": "no_hop", "product": None}], frm="cis")
    assert a["frac_switch"] is None and a["ci_low"] is None


def test_duplicate_traj_raises():                    # Review Focus 5: CSV antigo + rerun do mesmo TRAJ
    rows = [{"traj": "TRAJ7", "status": "died_S0", "product": None},
            {"traj": "TRAJ7", "status": "hop_ok", "product": "trans"}]
    try:
        aggregate(rows, frm="cis")
        assert False, "should raise"
    except ValueError as e:
        assert "duplicate traj: TRAJ7" in str(e)


def test_not_run_excluded():
    rows = [{"traj": "T1", "status": "hop_ok", "product": "trans"}, {"traj": "T2", "status": "not_run", "product": None}]
    a = aggregate(rows, frm="cis")
    assert a["n_total"] == 1 and a["n_hop_ok"] == 1


if __name__ == "__main__":
    fails = 0
    for n, f in sorted(globals().items()):
        if n.startswith("test_") and callable(f):
            try: f(); print("PASS", n)
            except Exception as e: fails += 1; print("FAIL", n, type(e).__name__, e)
    sys.exit(1 if fails else 0)
