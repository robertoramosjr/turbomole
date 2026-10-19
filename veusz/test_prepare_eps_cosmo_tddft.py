import prepare_eps_cosmo_tddft as m

def test_select_cover_and_cap():
    pop = {"a": .5, "b": .3, "c": .1, "d": .1}
    assert m.select(pop, set(pop), 0.8, 10) == ["a", "b"]
    assert m.select(pop, set(pop), 0.99, 3) == ["a", "b", "c"]

def test_cosmo_block(tmp_path):
    p = tmp_path / "coord"
    p.write_text("$coord\n 0 0 0 o\n 1 0 0 c\n 2 0 0 c\n 3 0 0 h\n$end\n")
    b = m.cosmo_block(str(p))
    assert "o  1-1" in b and "c  2-3" in b and "h  4-4" in b and b.startswith("$cosmo\n")
