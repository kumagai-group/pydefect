# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import numpy as np
import pytest
from pydefect.analyzer.defect_charge_info import DefectChargeInfo, ChargeDist
from pymatgen.electronic_structure.core import Spin
from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp import Chgcar
from vise.tests.helpers.assertion import assert_json_roundtrip


@pytest.fixture
def charge_dist():
    return ChargeDist(charge_center=(0.1, 0.1, 0.1), radial_dist=[0.5, 0.25, 0.25])


@pytest.fixture
def defect_charge_info(charge_dist):
    return DefectChargeInfo(distance_bins=[0.0, 1.0, 2.0, 2.5],
                            band_idxs=[10],
                            charge_dists=[[charge_dist]],
                            ave_charge_density=0.5)


def test_charge_dist_json_roundtrip(charge_dist, tmpdir):
    assert_json_roundtrip(charge_dist, tmpdir)


def test_defect_charge_info_json_roundtrip(defect_charge_info, tmpdir):
    assert_json_roundtrip(defect_charge_info, tmpdir)

#def test_bin_me(defect_charge_info, tmpdir):
#
# def test_make_defect_charge_info_is_charge_localized(defect_charge_info):
#     actual = defect_charge_info.is_charge_localized(10, Spin.up, threshold=(3 / 4 - 0.1))
#     assert actual is True
#     actual = defect_charge_info.is_charge_localized(10, Spin.up, threshold=(3 / 4 + 0.1))
#     assert actual is False
#
#
# def test_make_defect_charge_info_defect_radius(defect_charge_info):
#     assert defect_charge_info.defect_radius(10, Spin.up) == 0.5 * 0.5 + .25 * 1.5 + .25 * 2.5



"""
TODO
* DefectChargeInfo dataclass
- center
- radial_dist (-> charge radius)
- defect_region_radius

* make DefectChargeInfo

"""