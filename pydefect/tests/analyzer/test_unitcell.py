# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np
import pytest

from pydefect.analyzer.unitcell import Unitcell
from pydefect.tests.helpers.assertion import assert_msonable


@pytest.fixture
def unitcell():
    ele = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    ion = np.array([[10, 20, 30], [40, 50, 60], [70, 80, 90]])
    return Unitcell(vbm=0.1, cbm=5.1,
                    ele_dielectric_const=ele, ion_dielectric_const=ion)


def test_unitcell(unitcell):
    expected = np.array([[11, 22, 33], [44, 55, 66], [77, 88, 99]])
    np.testing.assert_array_almost_equal(unitcell.dielectric_constant, expected)


def test_msonable(unitcell):
    assert_msonable(unitcell)

