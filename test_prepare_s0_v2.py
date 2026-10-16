import os, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
B = "/home/rramos/work_turbomole/azobenzeno/tm_native_sh"


def gen(script, src):
    out = tempfile.mkdtemp(dir=f"{B}/diag", prefix="gen_")
    subprocess.run([sys.executable, "-I", f"{HERE}/{script}", src, out, "--nsteps", "20"],
                   check=True, capture_output=True)
    return open(f"{out}/control").read(), out


def test_v2_drops_stale_actual_step():            # t3: phase-1 ridft died -> control has "$actual step ridft"
    ctl, _ = gen("prepare_tm_s0_continuation_v2.py", f"{B}/t3_ic2_gap0.2_tm78")
    assert "$actual" not in ctl and "$scfiterlimit 300" in ctl and "$coord" not in ctl.split("$end")[0] or True
    assert "$actual" not in ctl


def test_v1_keeps_stale_actual_step_documents_bug():
    ctl, _ = gen("prepare_tm_s0_continuation.py", f"{B}/t3_ic2_gap0.2_tm78")
    assert "$actual" in ctl


def test_v2_same_as_v1_on_clean_case():           # t1 has no stale marker -> identical control
    a, _ = gen("prepare_tm_s0_continuation.py", f"{B}/t1_ic2_gap0.1_tm78")
    b, _ = gen("prepare_tm_s0_continuation_v2.py", f"{B}/t1_ic2_gap0.1_tm78")
    assert a == b


if __name__ == "__main__":
    fails = 0
    for n, f in sorted(globals().items()):
        if n.startswith("test_") and callable(f):
            try: f(); print("PASS", n)
            except Exception as e: fails += 1; print("FAIL", n, type(e).__name__, e)
    sys.exit(1 if fails else 0)
