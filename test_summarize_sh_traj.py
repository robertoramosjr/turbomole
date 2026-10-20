import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarize_sh_traj import summarize

B = "/home/rramos/work_turbomole/azobenzeno/tm_native_sh"


def test_cis_traj1_hop_trans():
    r = summarize(f"{B}/batch_cis420_gap0.2/TRAJ1")
    assert r["status"] == "hop_ok" and r["product"] == "trans"
    assert 20 <= r["t_hop_fs"] <= 40 and 75 <= abs(r["cnnc_hop"]) <= 95


def test_trans_traj2_died_s1():
    r = summarize(f"{B}/batch_trans425_gap0.2/TRAJ2")
    assert r["status"] == "died_S1" and r["product"] is None


def test_trans_traj1_hop():
    r = summarize(f"{B}/batch_trans425_gap0.2/TRAJ1")
    assert r["status"] == "hop_ok" and r["product"] == "trans"


def test_t3_died_s0():
    r = summarize(f"{B}/t3_ic2_gap0.2_tm78", s0_dir=f"{B}/s0_cont/t3_ic2_gap0.2_tm78")
    assert r["status"] == "died_S0" and r["product"] is None
    assert r["n_hops"] == 3


def test_missing_dir_raises():                       # P05: caminho errado nao pode virar no_hop
    try:
        summarize(f"{B}/nao_existe/TRAJ1")
        assert False, "should raise"
    except FileNotFoundError:
        pass


def test_virgin_dir_not_run():                       # TRAJ7 do lote v2 nunca rodou (lote cancelado)
    r = summarize(f"{B}/batch_trans425v2_gap0.2/TRAJ7")
    assert r["status"] == "not_run" and r["product"] is None


if __name__ == "__main__":  # sem pytest na access: executor mínimo
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn(); print("PASS", name)
            except Exception as e:
                fails += 1; print("FAIL", name, type(e).__name__, e)
    sys.exit(1 if fails else 0)
