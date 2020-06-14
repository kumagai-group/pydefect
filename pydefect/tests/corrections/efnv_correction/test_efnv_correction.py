# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.corrections.efnv_correction.efnv_correction import \
    ExtendedFnvCorrection, PotentialSite
from pydefect.tests.helpers.assertion import assert_msonable


@pytest.fixture
def efnv_correction():
    s1 = PotentialSite(specie="H", distance=1.999, potential=1.0, pc_potential=None)
    s2 = PotentialSite(specie="He", distance=2.0001, potential=1.5, pc_potential=0.2)
    s3 = PotentialSite(specie="He", distance=3.0, potential=2.0, pc_potential=0.3)

    return ExtendedFnvCorrection(charge=10,
                                 point_charge_correction=1.0,
                                 defect_region_radius=2.0,
                                 sites=[s1, s2, s3],
                                 additional_correction=100)


def test_extended_fnv_correction_msonable(efnv_correction):
    assert_msonable(efnv_correction)


def test_extended_fnv_correction_ave_pot_diff(efnv_correction):
    assert efnv_correction.average_potential_diff == ((1.5 - 0.2) + (2.0 - 0.3)) / 2


def test_extended_fnv_correction_alignment_correction(efnv_correction):
    assert efnv_correction.alignment_correction == -10 * ((1.5 - 0.2) + (2.0 - 0.3)) / 2


def test_extended_fnv_correction_correction_energy(efnv_correction):
    assert efnv_correction.correction_energy == 1.0 - 10 * ((1.5 - 0.2) + (2.0 - 0.3)) / 2 + 100


def test_defect_site_diff_pot():
    s = PotentialSite(specie="H", distance=1.999, potential=1.0, pc_potential=0.1)
    assert s.diff_pot == 1.0 - 0.1


def test_defect_site_msonable(efnv_correction):
    s = PotentialSite(specie="H", distance=1.999, potential=1.0, pc_potential=0.1)
    assert_msonable(s)

