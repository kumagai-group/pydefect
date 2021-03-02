# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from numpy import pi, sqrt, exp, cos
from pydefect.corrections.ewald import Ewald
from scipy.special import erfc
from vise.tests.helpers.assertion import assert_msonable


@pytest.fixture
def ewald():
    return Ewald(lattice=np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]]),
                 dielectric_tensor=np.array([[3, 0, 0], [0, 3, 0], [0, 0, 3]]),
                 ewald_param=0.1,
                 accuracy=0.1)


def test_msonable(ewald):
    assert_msonable(ewald)


def test_ewald(ewald):
    np.testing.assert_almost_equal(
        ewald.lattice, np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]]))
    np.testing.assert_almost_equal( # 0.5 * 2 * pi
        ewald.rec_lattice, np.array([[pi, 0, 0], [0, pi, 0], [0, 0, pi]]))
    np.testing.assert_almost_equal(ewald.volume, 8.0)
    np.testing.assert_almost_equal(ewald.cube_root_vol, 2.0)
    np.testing.assert_almost_equal(ewald.det_epsilon, 27)
    np.testing.assert_almost_equal(ewald.root_epsilon, 3 * np.sqrt(3))
    np.testing.assert_almost_equal(
        ewald.epsilon_inv, np.array([[1/3, 0, 0], [0, 1/3, 0], [0, 0, 1/3]]))
    np.testing.assert_almost_equal(
        ewald.mod_ewald_param, 0.1 / 2.0 * 3 * np.sqrt(3))
    np.testing.assert_almost_equal(
        ewald.diff_pot, -0.25 / 8.0 / (0.1 / 2.0 * 3 * np.sqrt(3)) ** 2)


def test_ewald_vector_nums():
    ewald = Ewald(lattice=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                  dielectric_tensor=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                  ewald_param=1,
                  accuracy=1.0001)
    assert ewald.r_vector_nums == [2, 2, 2]

    ewald = Ewald(lattice=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                  dielectric_tensor=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                  ewald_param=1,
                  accuracy=6.28318531 / 2 + 0.0001)
    assert ewald.g_vector_nums == [2, 2, 2]


def test_ewald_xyz(ewald):
    actual = ewald.xyz([1, 1, 1], shift=[0.1, 0.1, 0.1])
    assert len(actual) == 27
    np.testing.assert_almost_equal(actual[0], np.array([-1.1, -1.1, -1.1]))


def test_ewald_lattice_set(ewald):
    np.testing.assert_almost_equal(ewald.r_lattice_set()[0], np.array([-2, -2, -2]))
    np.testing.assert_almost_equal(ewald.g_lattice_set()[0], np.array([-3.14159265, -3.14159265, -3.14159265]))


def test_ewald_real(ewald, mocker):
    mock = mocker.patch("pydefect.corrections.ewald.Ewald.r_lattice_set")
    mock.return_value = np.array([[-2, -2, -2]])
    root_r_inv_epsilon_r = sqrt(2 / 3 * 2 * 3)
    actual = ewald.ewald_real(include_self=False, shift=[0, 0, 0])
    expected = erfc(ewald.mod_ewald_param * root_r_inv_epsilon_r) / root_r_inv_epsilon_r / (4 * pi * ewald.root_epsilon)
    assert actual == expected


def test_ewald_rec(ewald, mocker):
    mock = mocker.patch("pydefect.corrections.ewald.Ewald.g_lattice_set")
    mock.return_value = np.array([[2, 2, 2]])
    g_epsilon_g = 2 * 3 * 2 * 3
    actual = ewald.ewald_rec([0.2, 0.2, 0.2])
    cos_term = cos(2 * 2 * 0.2 * 3)
    expected = exp(-g_epsilon_g / 4 / ewald.mod_ewald_param ** 2) / g_epsilon_g * cos_term / ewald.volume
    assert actual == expected

#
# def test_ewald_speed():
#     ewald = Ewald(lattice=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
#                   dielectric_tensor=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
#                   accuracy=30)
#     assert ewald.ewald_real(include_self=False, shift=None) == 0.11400074023354886



"""
TODO
- Create a container Ewald class
- 

DONE
"""