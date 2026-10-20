import os, re, shutil, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
B = "/home/rramos/work_turbomole/azobenzeno/tm_native_sh"
GROUPS = ("population_el", "prob_fssh", "sh_coeffs", "egradmon")


def seed(src, maxlog):
    d = tempfile.mkdtemp(dir=f"{B}/diag", prefix="seed_")
    ctl = f"{d}/control"
    shutil.copy(src, ctl)
    subprocess.run([sys.executable, "-I", f"{HERE}/seed_sh_log_groups.py", ctl, "--maxlog", str(maxlog)],
                   check=True, capture_output=True)
    return open(ctl).read(), ctl, d


def test_seeds_four_groups():                      # TRAJ8 do lote v2 e' virgem: sem os grupos
    src = f"{B}/batch_trans425v2_gap0.2/TRAJ8/control"
    assert "$population_el" not in open(src).read()
    ctl, _, _ = seed(src, 20)
    assert len(re.findall(r"^ maxlog= *20$", ctl, re.M)) == 4
    for g in GROUPS:
        assert re.search(rf"^\${g} log=true$", ctl, re.M)
        assert f" logbase={g}log" in ctl


def test_idempotent():
    src = f"{B}/batch_trans425v2_gap0.2/TRAJ8/control"
    ctl1, path, d = seed(src, 1000)
    subprocess.run([sys.executable, "-I", f"{HERE}/seed_sh_log_groups.py", path, "--maxlog", "1000"], check=True)
    assert open(path).read() == ctl1


def test_replaces_existing():                       # diag/maxlog20/control ja tem os grupos com 256
    src = f"{B}/diag/maxlog20/control"
    n0 = len(open(src).read().splitlines())
    ctl, _, _ = seed(src, 1000)
    assert not re.search(r"maxlog= *256\b", ctl.split("$title")[0])
    assert len(re.findall(r"^ maxlog= *1000$", ctl, re.M)) == 4
    assert len(ctl.splitlines()) == n0


if __name__ == "__main__":
    fails = 0
    for n, f in sorted(globals().items()):
        if n.startswith("test_") and callable(f):
            try: f(); print("PASS", n)
            except Exception as e: fails += 1; print("FAIL", n, type(e).__name__, e)
    sys.exit(1 if fails else 0)
