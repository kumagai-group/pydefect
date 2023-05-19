# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

import numpy as np
import pytest

from pydefect.analyzer.unitcell import Unitcell
from vise.tests.helpers.assertion import assert_yaml_roundtrip


@pytest.fixture
def unitcell():
    ele = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [0, 8.0, 9.0]]
    ion = [[10.0, 20.0, 30.0], [40.0, 50.0, 60.0], [0, 80.0, 90.0]]
    e_mass = [[11.0, 12.0, 13.0], [24.0, 25.0, 26.0], [37.0, 38.0, 39.0]]
    h_mass = [[111.0, 112.0, 113.0], [124.0, 125.0, 126.0], [137.0, 138.0, 139.]]
    return Unitcell(system="MgAl2O4", vbm=-1.5053, cbm=5.1,
                    ele_dielectric_const=ele, ion_dielectric_const=ion,
                    electron_mass=e_mass, hole_mass=h_mass)


def test_unitcell(unitcell):
    expected = np.array([[11, 22, 33], [44, 55, 66], [0, 88, 99]])
    np.testing.assert_array_almost_equal(unitcell.dielectric_constant, expected)


def test_average_diele_const(unitcell):
    assert unitcell.ave_ele_diele == (1.0 + 5.0 + 9.0) / 3
    assert unitcell.ave_diele == (11.0 + 55.0 + 99.0) / 3


def test_average_masses(unitcell):
    assert unitcell.ave_ele_mass == (11.0+25.0+39.0) / 3
    assert unitcell.ave_hole_mass == (111.0+125.0+139.0) / 3


def test_effective_ionic_diele_const(unitcell):
    assert unitcell.effective_ionic_diele_const[0][0] == 1/(1/1.0 - 1/11.0)
    assert unitcell.effective_ionic_diele_const[2][0] == 0.0


def test_yaml_roundtrip(unitcell, tmpdir):
    expected_text = """system: MgAl2O4
vbm: -1.5053
cbm: 5.1
ele_dielectric_const:
- - 1.0
  - 2.0
  - 3.0
- - 4.0
  - 5.0
  - 6.0
- - 0
  - 8.0
  - 9.0
ion_dielectric_const:
- - 10.0
  - 20.0
  - 30.0
- - 40.0
  - 50.0
  - 60.0
- - 0
  - 80.0
  - 90.0
electron_mass:
- - 11.0
  - 12.0
  - 13.0
- - 24.0
  - 25.0
  - 26.0
- - 37.0
  - 38.0
  - 39.0
hole_mass:
- - 111.0
  - 112.0
  - 113.0
- - 124.0
  - 125.0
  - 126.0
- - 137.0
  - 138.0
  - 139.0
"""
    assert_yaml_roundtrip(unitcell, tmpdir, expected_text)
