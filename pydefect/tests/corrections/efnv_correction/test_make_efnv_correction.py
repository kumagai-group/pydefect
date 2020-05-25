# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np

from pydefect.corrections.efnv_correction.make_efnv_correction import \
    calc_max_sphere_radius


# @pytest.fixture
# def pot():
    # ewald = Ewald(lattice=np.array([[10, 0, 0], [0, 10, 0], [0, 0, 10]]),
    #               dielectric_tensor=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    #               accuracy=5)
    # return AnisotropicPotential(defect_coord=np.array([0.0, 0.0, 0.0]),
    #                             target_coords=np.array([[0.5, 0.5, 0.5]]),
    #                             ewald=ewald)


# def test_pot(pot):
#     assert pot.lattice_energy * 180.95128169876497 == -2.0428038899554415



def test_calc_max_sphere_radius():
    lattice_vectors_1 = np.array([[5, 0, 0], [0, 10, 0], [0, 0, 20]])
    lattice_vectors_2 = np.array([[10, 0, 0], [0, 10, 0], [10, 10, 10]])
    radius_1 = calc_max_sphere_radius(lattice_vectors_1)
    radius_2 = calc_max_sphere_radius(lattice_vectors_2)
    assert radius_1 == 10.0
    assert radius_2 == 5.0

