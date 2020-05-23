# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from pymatgen import Lattice

from pydefect.corrections.efnv_correction.efnv_correction import \
    create_lattice_set
from pydefect.corrections.efnv_correction.ewald import Ewald, \
    calc_max_sphere_radius, create_lattice_set_2
from pydefect.tests.helpers.assertion import assert_msonable


def test_calc_max_sphere_radius():
    lattice_vectors_1 = np.array([[5, 0, 0], [0, 10, 0], [0, 0, 20]])
    assert calc_max_sphere_radius(lattice_vectors_1) == 10.0

    lattice_vectors_2 = np.array([[10, 0, 0], [0, 10, 0], [10, 10, 10]])
    assert calc_max_sphere_radius(lattice_vectors_2) == 5.0


lattice_vectors_1 = np.array([[10, 10, 0], [10, -10, 0], [0, 0, 14]])
def test_lattice_set():
    actual = create_lattice_set(lattice_vectors_1, max_length=15)
    expected = [[-10, -10, 0], [-10, 10, 0], [0, 0, -14], [0, 0, 0], [0, 0, 14],
                [10, -10, 0], [10, 10, 0]]
    assert actual == expected


def test_lattice_speed_set():
    actual = create_lattice_set_2(lattice_vectors_1, max_length=1500)
    print(len(actual))


@pytest.fixture
def ewald():
    return Ewald(lattice=Lattice.cubic(3.0).matrix,
                 dielectric_tensor=np.array([[1.0, 0.0, 0.0],
                                             [0.0, 1.0, 0.0],
                                             [0.0, 0.0, 1.0]]),
                 ewald_param=1.0,
                 accuracy=1.0)


def test_ewald(ewald):
    assert ewald.volume == 3.0 ** 3
    assert (ewald.reciprocal_lattice == np.eye(3) * (2 * np.pi / 3)).all()


def test_msonable(ewald):
    assert_msonable(ewald)

# def test_make_ewald(ewald):



"""
TODO
- Create a container Ewald class
- 

DONE
"""