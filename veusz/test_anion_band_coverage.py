import numpy as np
import anion_band_coverage as m

def test_lambda_max_single_stick(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "sticks", lambda p: [(4.0, 1.0)])
    grid = np.linspace(2.0, 6.0, 4000)
    lam, tot, emax = m.lambda_max({"j"}, {"j": 0.5}, "x", 0.4, grid)
    assert abs(lam - m.HC / 4.0) < 0.2 and tot == 0.5
