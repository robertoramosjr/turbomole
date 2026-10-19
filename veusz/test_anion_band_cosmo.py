import os
import numpy as np
import anion_band_cosmo as m

ROOT = os.path.expanduser("~/work_turbomole/artepillin_C_studies")


def test_single_stick_lambda(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "sticks", lambda p: [(4.0, 0.5)])
    assert abs(m.lambda_max(["x"], [1.0], 0.2, 3.0, 4.96) - m.HC / 4.0) < 0.2


def test_window_ignores_far_root(monkeypatch):
    monkeypatch.setattr(m, "sticks", lambda p: [(4.0, 0.5), (6.9, 0.9)])  # 6.9 eV ~ 178 nm artefact
    assert abs(m.lambda_max(["x"], [1.0], 0.2, 3.0, 4.96) - m.HC / 4.0) < 0.2


def test_run_has_both_media():
    rows = m.run(ROOT, 0.4, 3.0, 4.96)
    assert {(r["state"], r["medium"]) for r in rows} == {("N", "vacuum"), ("N", "cosmo"), ("A", "vacuum"), ("A", "cosmo")}
    assert all(r["n_conf"] > 0 for r in rows)
