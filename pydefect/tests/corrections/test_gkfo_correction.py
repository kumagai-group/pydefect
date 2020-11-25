# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from pydefect.corrections.abstract_correction import Correction
from pydefect.corrections.efnv_correction import PotentialSite
from pydefect.corrections.gkfo_correction import GkfoCorrection


@pytest.fixture
def gkfo_correction(efnv_correction):
    s1 = PotentialSite(
        specie="H", distance=1.999, potential=10, pc_potential=None)
    s2 = PotentialSite(
        specie="He", distance=2.0001, potential=15, pc_potential=2)
    s3 = PotentialSite(
        specie="He", distance=3.0, potential=20, pc_potential=3)

    return GkfoCorrection(init_efnv_correction=efnv_correction,
                          additional_charge=1,
                          pc_2nd_term=10.0,
                          gkfo_sites=[s1, s2, s3],
                          ave_dielectric_tensor=4.0,
                          ave_electronic_dielectric_tensor=2.0)


def test_is_correction(gkfo_correction):
    assert isinstance(gkfo_correction, Correction)


def test_repr(gkfo_correction):
    expected = """+--------------------+-------+
| charge             |  10   |
| additional charge  |   1   |
| pc 1st term        |   0.2 |
| pc 2nd term        |  10   |
| alignment 1st term | -15   |
| alignment 2nd term |  -1.5 |
| alignment 3rd term | -75   |
| correction energy  | -81.3 |
+--------------------+-------+"""
    assert repr(gkfo_correction) == expected


def test_gkfo_correction_charge(gkfo_correction):
    assert gkfo_correction.charge == 10


def test_gkfo_correction_defect_region_radius(gkfo_correction):
    assert gkfo_correction.defect_region_radius == 2.0


def test_gkfo_correction_pc_1st_term(gkfo_correction):
    assert gkfo_correction.pc_1st_term == 1.0 * 2 / 10 * 1


def test_gkfo_correction_sum_pc_correction(gkfo_correction):
    assert gkfo_correction.sum_pc_correction == 1.0 * 2 / 10 * 1 + 10.0


def test_extended_fnv_correction_alignment_correction(gkfo_correction):
    assert gkfo_correction.average_potential_diff_by_addition == ((15 - 2) + (20 - 3)) / 2
    assert gkfo_correction.alignment_1st_term == - ((15 - 2) + (20 - 3)) / 2
    assert gkfo_correction.alignment_2nd_term == - ((1.5 - 0.2) + (2 - 0.3)) / 2
    assert gkfo_correction.alignment_3rd_term == - ((15 - 2) + (20 - 3)) / 2 * 10 * 2 / 4
    assert gkfo_correction.sum_alignment_term == - 15 - 1.5 - 75


def test_extended_fnv_correction_correction_energy(gkfo_correction):
    assert gkfo_correction.correction_energy == 0.2 + 10.0 - 15 - 1.5 - 75


