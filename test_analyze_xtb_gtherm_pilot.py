import os, sys, tempfile
sys.path.insert(0, os.path.dirname(__file__))
import analyze_xtb_gtherm_pilot as a


def test_parse_g_rrho():
    d = tempfile.mkdtemp()
    p = os.path.join(d, "bhess.out")
    open(p, "w").write("   :: G(RRHO) contrib.           0.343116523750 Eh   ::\n")
    assert abs(a.parse_g_rrho(p) - 0.343116523750) < 1e-12
    assert a.parse_g_rrho(os.path.join(d, "none.out")) is None


def test_compare_perfect_and_offset():
    h = [0.340, 0.341, 0.3425]
    x = [v + 0.1 for v in h]            # constant offset must cancel
    m = a.compare(h, x)
    assert m["max_abs"] < 1e-9 and abs(m["pearson"] - 1) < 1e-9


def test_compare_detects_error():
    h = [0.340, 0.341, 0.342]
    x = [0.340, 0.3415, 0.342]          # 0.0005 Eh = 0.31 kcal/mol error on the middle one
    assert abs(a.compare(h, x)["max_abs"] - 0.0005 * a.HARTREE_TO_KCAL) < 1e-6
