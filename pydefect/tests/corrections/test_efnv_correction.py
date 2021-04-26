# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pydefect.corrections.abstract_correction import Correction

from pydefect.corrections.efnv_correction import \
    PotentialSite
from vise.tests.helpers.assertion import assert_msonable


def test_extended_fnv_correction_msonable(efnv_correction):
    assert isinstance(efnv_correction, Correction)


def test_repr(efnv_correction):
    expected = """+-------------------+-----+
| charge            |  10 |
| pc term           |   1 |
| alignment term    | -15 |
| correction energy | -14 |
+-------------------+-----+"""
    assert str(efnv_correction) == expected


def test_extended_fnv_correction_ave_pot_diff(efnv_correction):
    expected = ((1.5 - 0.2) + (2.0 - 0.3)) / 2
    assert efnv_correction.average_potential_diff == expected


def test_extended_fnv_correction_alignment_correction(efnv_correction):
    expected = -10 * ((1.5 - 0.2) + (2.0 - 0.3)) / 2
    assert efnv_correction.alignment_correction == expected


def test_extended_fnv_correction_correction_energy(efnv_correction):
    expected = 1.0 - 10 * ((1.5 - 0.2) + (2.0 - 0.3)) / 2
    assert efnv_correction.correction_energy == expected


def test_extended_fnv_correction_correction_dict(efnv_correction):
    expected = {"pc term": 1.0, "alignment term": -15.0}
    assert efnv_correction.correction_dict == expected


def test_defect_site_diff_pot():
    s = PotentialSite("H", distance=1.999, potential=1.0, pc_potential=0.1)
    assert s.diff_pot == 1.0 - 0.1


def test_defect_site_msonable(efnv_correction):
    assert_msonable(PotentialSite(specie="H", distance=1.999, potential=1.0,
                                  pc_potential=0.1))

