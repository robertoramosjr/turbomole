import math
import pka_estimate as m


def test_macro_matches_closure():
    p = {"N_to_A": 4.0, "N_to_B": 8.0, "A_to_D": 11.0, "B_to_D": 7.0}
    pka1, pka2 = m.macro(p)
    assert abs(pka1 - 4.0) < 1e-3                      # B negligible
    assert abs(pka1 + pka2 - 15.0) < 1e-3              # Ka1*Ka2 = KNA*KAD


def test_analogue_additive():
    refs = [{"name": "26-dimethylphenol", "exp_pka": "10.61"}, {"name": "phenol", "exp_pka": "9.99"},
            {"name": "4-coumarate_phenol", "exp_pka": "9.45"}]
    assert abs(m.analogue(refs) - 10.07) < 1e-9
    assert abs(m.analogue(refs, 9.92) - 10.54) < 1e-9


def test_quad():
    assert abs(m.quad(3, 4) - 5) < 1e-12
