# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from pymatgen import Lattice

from pydefect.corrections.efnv_correction.ewald import Ewald
from pydefect.tests.helpers.assertion import assert_msonable


@pytest.fixture
def ewald():
    return Ewald(lattice=Lattice.cubic(3.0),
                 dielectric_tensor=np.array([[1.0, 0.0, 0.0],
                                             [0.0, 1.0, 0.0],
                                             [0.0, 0.0, 1.0]]),
                 ewald_param=1.0,
                 accuracy=1.0)


def test_ewald(ewald):
    assert ewald.volume == 3.0 ** 3
    assert (ewald.rec_lattice_matrix == np.eye(3) * (2 * np.pi / 3)).all()


def test_msonable(ewald):
    assert_msonable(ewald)


"""
TODO
- Create a container Ewald class
- 

DONE
"""