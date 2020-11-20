# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import numpy as np
from pydefect.analyzer.defect_charge_distribution import RadialDist
from pymatgen import Structure, Lattice, Spin
from pymatgen.io.vasp import Chgcar


def test():
    structure = Structure(Lattice.cubic(3), ["H"], [[0, 0, 0]])
    data = {"total": np.array([[[3]*3]*3]*3), "diff": np.array([[[-1]*3]*3]*3)}
    chgcar = Chgcar(structure, data)

    rad = RadialDist(chgcar, [0, 0, 0])
    assert len(rad.distances_data) == 1 + 6 + 12
    np.testing.assert_almost_equal(rad.distances_data[0][0], np.array([-1/3, -1./3, 0.0]))
    assert rad.distances_data[0][1] == np.sqrt(2)
    assert rad.distances_data[0][2] == 24
    np.testing.assert_almost_equal(rad.distances_data[0][3], np.array([-1, -1, 0.0]))

    hist_data, half_point, summed = rad.histogram(Spin.up)
    np.testing.assert_almost_equal(hist_data[0], [5.00000000e-02, 1.16355283e-03])


