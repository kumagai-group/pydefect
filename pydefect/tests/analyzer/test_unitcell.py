# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

import numpy as np
import pytest

from pydefect.analyzer.unitcell import Unitcell
from vise.tests.helpers.assertion import assert_yaml_roundtrip


@pytest.fixture
def unitcell():
    ele = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
    ion = [[10.0, 20.0, 30.0], [40.0, 50.0, 60.0], [70.0, 80.0, 90.0]]
    return Unitcell(system="MgAl2O4", vbm=-1.5053, cbm=5.1,
                    ele_dielectric_const=ele, ion_dielectric_const=ion)


def test_unitcell(unitcell):
    expected = np.array([[11, 22, 33], [44, 55, 66], [77, 88, 99]])
    np.testing.assert_array_almost_equal(unitcell.dielectric_constant, expected)


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
- - 7.0
  - 8.0
  - 9.0
ion_dielectric_const:
- - 10.0
  - 20.0
  - 30.0
- - 40.0
  - 50.0
  - 60.0
- - 70.0
  - 80.0
  - 90.0
"""
    assert_yaml_roundtrip(unitcell, tmpdir, expected_text)
