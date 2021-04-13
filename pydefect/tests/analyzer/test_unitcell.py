# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

import numpy as np
import pytest

from pydefect.analyzer.unitcell import Unitcell
from vise.tests.helpers.assertion import assert_msonable, assert_json_roundtrip


@pytest.fixture
def unitcell():
    ele = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
    ion = [[10.0, 20.0, 30.0], [40.0, 50.0, 60.0], [70.0, 80.0, 90.0]]
    return Unitcell(vbm=0.1, cbm=5.1,
                    ele_dielectric_const=ele, ion_dielectric_const=ion)


def test_unitcell(unitcell):
    expected = np.array([[11, 22, 33], [44, 55, 66], [77, 88, 99]])
    np.testing.assert_array_almost_equal(unitcell.dielectric_constant, expected)


def test_json_round_trip(unitcell, tmpdir):
    assert_json_roundtrip(unitcell, tmpdir)


def test_yaml_roundtrip(unitcell, tmpdir):
    tmpdir.chdir()
    unitcell.to_yaml("a.yaml")
    expected = """vbm: 0.1
cbm: 5.1
ele_dielectric_const:
- [1.0, 2.0, 3.0]
- [4.0, 5.0, 6.0]
- [7.0, 8.0, 9.0]
ion_dielectric_const:
- [10.0, 20.0, 30.0]
- [40.0, 50.0, 60.0]
- [70.0, 80.0, 90.0]
"""
    assert Path("a.yaml").read_text() == expected

    actual = Unitcell.from_yaml("a.yaml").as_dict()
    expected = unitcell.as_dict()
    for k, v in actual.items():
        assert v == expected[k]
