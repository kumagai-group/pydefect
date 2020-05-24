# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest

from pydefect.corrections.efnv_correction.ewald import Ewald
from pydefect.corrections.efnv_correction.potential import AnisotropicPotential


@pytest.fixture
def pot():
    ewald = Ewald(lattice=np.array([[10, 0, 0], [0, 10, 0], [0, 0, 10]]),
                  dielectric_tensor=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                  accuracy=5)
    return AnisotropicPotential(defect_coord=np.array([0.0, 0.0, 0.0]),
                                target_coords=np.array([[0.5, 0.5, 0.5]]),
                                ewald=ewald)


def test_pot(pot):
    assert pot.lattice_energy * 180.95128169876497 == -2.0428038899554415
