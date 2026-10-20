import filecmp, os, re, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
B = "/home/rramos/work_turbomole/azobenzeno/tm_native_sh"
SRC = f"{B}/batch_cis420_gap0.2/TRAJ1"
FS = 41.341373


def gen(script, *extra):
    out = tempfile.mkdtemp(dir=f"{B}/diag", prefix="gen3_")
    r = subprocess.run([sys.executable, "-I", f"{HERE}/{script}", SRC, out, "--nsteps", "20", *extra],
                       check=True, capture_output=True, text=True)
    return out, r.stdout


def header_t(out):
    return float(re.search(r"^t=\s*(\S+)", open(f"{out}/mdlog.1").read(), re.M).group(1))


def test_v3_k0_equals_v2():
    a, _ = gen("prepare_tm_s0_continuation_v2.py")
    b, _ = gen("prepare_tm_s0_continuation_v3.py")
    for f in ("control", "coord", "mdlog.1", "mdmaster"):
        assert filecmp.cmp(f"{a}/{f}", f"{b}/{f}", shallow=False), f


def test_frame_back_2_time_and_block():
    a, _ = gen("prepare_tm_s0_continuation_v3.py")
    b, _ = gen("prepare_tm_s0_continuation_v3.py", "--frame-back", "2")
    assert abs((header_t(a) - header_t(b)) - 2 * 0.5 * FS) < 1e-3
    body = open(f"{b}/mdlog.1").read().split("$current\n")[1].split("$end")[0].strip("\n").split("\n")
    assert len(body) == 1 + 2 * 24 + 1                       # t= + 24 coords + dt + 24 velocities
    assert body[0].startswith("t=")


def test_mos_from():
    alt = f"{B}/diag/s0_matrix_mos_probe"
    open(alt, "w").write("$scfmo marker\n$end\n")
    b, _ = gen("prepare_tm_s0_continuation_v3.py", "--mos-from", alt)
    assert filecmp.cmp(f"{b}/mos", alt, shallow=False)


if __name__ == "__main__":
    fails = 0
    for n, f in sorted(globals().items()):
        if n.startswith("test_") and callable(f):
            try: f(); print("PASS", n)
            except Exception as e: fails += 1; print("FAIL", n, type(e).__name__, e)
    sys.exit(1 if fails else 0)
