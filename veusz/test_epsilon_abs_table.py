import numpy as np
from epsilon_abs_table import eps_from_profile, fwhm_ev
from plot_paper2_figures import HC


def test_single_stick():
    g = np.linspace(2, 7, 5001)
    prof = np.exp(-((g - 4.0) ** 2) / (2 * 0.4 ** 2))  # f = 1 at 4 eV
    assert abs(eps_from_profile(prof, g, 0.4, (3, 5.5)) - 4.05e4) < 0.01e4


def test_fwhm_gaussian():
    e = np.linspace(3, 5, 4001)
    wl = HC / e
    y = np.exp(-((e - 4) ** 2) / (2 * 0.3 ** 2))
    fw, note = fwhm_ev(wl, y)
    assert abs(fw - 0.3 * 2.3548) < 0.01 and note == ""


def test_subset_spectrum_weights(monkeypatch):
    import epsilon_abs_table as m
    monkeypatch.setattr(m, "sticks", lambda p: [(4.0, 1.0)])
    g = np.linspace(2, 7, 5001)
    prof = m.subset_spectrum(["a", "b"], [0.2, 0.6], g, 0.4)  # weights renormalised -> sum f = 1
    assert abs(m.eps_from_profile(prof, g, 0.4, (3, 5.5)) - 4.05e4) < 0.01e4
