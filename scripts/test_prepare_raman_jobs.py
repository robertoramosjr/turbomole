import prepare_raman_jobs as m


def test_strip_tddft_keeps_rest():
    ctl = ["$dft", "    functional   pbe0", "$rik", "$scfinstab rpas", "$soes", " a           10",
           "$rpacor 500", "$xctype hybrid-GGA", "$end"]
    assert m.strip_tddft(ctl) == ["$dft", "    functional   pbe0", "$rik", "$xctype hybrid-GGA", "$end"]


def test_select_cover():
    pop = {"a": 0.6, "b": 0.3, "c": 0.1}
    sel, acc = m.select(pop, {"a", "b", "c"}, 0.85, 99)
    assert sel == ["a", "b"] and abs(acc - 0.9) < 1e-9
    assert m.select(pop, {"a", "b", "c"}, 0.99, 1)[0] == ["a"]
