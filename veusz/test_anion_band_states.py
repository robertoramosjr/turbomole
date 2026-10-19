import os
import anion_band_states as m

P = os.path.expanduser("~/work_turbomole/artepillin_C_studies/artepillin_c_monoanion_a_d3_hse06/layer_b/labeled_jobs/job_0062/escf.out")

def test_first_root():
    r = m.parse_roots(P)
    assert len(r) == 60 and r[0]["root"] == 1
    assert abs(r[0]["e"] - 2.5547) < 1e-3 and r[0]["f"] < 1e-3
    assert r[0]["c"][0][:3] == (81, -1.22, 82) and r[0]["c"][0][4] == 76.8
