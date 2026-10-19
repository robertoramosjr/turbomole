import numpy as np
import anion_band_geometry as m

def test_dihedral_and_twist():
    p = [np.array(v, float) for v in ([1, 1, 0], [0, 1, 0], [0, 0, 0], [0, 0, 1])]
    assert abs(abs(m.dihedral(*p)) - 90) < 1e-6
    assert abs(m.twist(170) - 10) < 1e-9 and m.twist(-5) == 5
