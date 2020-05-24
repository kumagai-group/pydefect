# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.corrections.efnv_correction.efnv_correction import \
    ExtendedFnvCorrection
from pydefect.tests.helpers.assertion import assert_msonable


@pytest.fixture
def efnv_correction():
    return ExtendedFnvCorrection(charge=10,
                                 point_charge_correction=1.0,
                                 defect_region_radius=2.0,
                                 species=["H", "He", "He"],
                                 distances=[1.999, 2.0001, 3.0],
                                 site_ele_pot=[1.0, 1.5, 2.0],
                                 point_charge_pot=[0.1, 0.2, 0.3],
                                 length_unit="Å",
                                 potential_unit="V",
                                 additional_correction=100)


def test_extended_fnv_correction_msonable(efnv_correction):
    assert_msonable(efnv_correction)


def test_extended_fnv_correction_diff_pot(efnv_correction):
    assert efnv_correction.diff_pot == [1.0 - 0.1, 1.5 - 0.2, 2.0 - 0.3]


def test_extended_fnv_correction_ave_pot_diff(efnv_correction):
    assert efnv_correction.ave_pot_diff == ((1.5 - 0.2) + (2.0 - 0.3)) / 2


def test_extended_fnv_correction_alignment_correction(efnv_correction):
    assert efnv_correction.alignment_correction == -10 * ((1.5 - 0.2) + (2.0 - 0.3)) / 2


def test_extended_fnv_correction_correction_energy(efnv_correction):
    assert efnv_correction.correction_energy == 1.0 - 10 * ((1.5 - 0.2) + (2.0 - 0.3)) / 2 + 100


"""
TODO
-

DONE
"""