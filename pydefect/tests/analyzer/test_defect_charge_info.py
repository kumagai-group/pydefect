# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import numpy as np
import pytest
from pydefect.analyzer.defect_charge_info import DefectChargeInfo, AveChargeDensityDist
from pymatgen.electronic_structure.core import Spin
from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp import Chgcar
from vise.tests.helpers.assertion import assert_json_roundtrip


@pytest.fixture
def ave_charge_density_dist():
    return AveChargeDensityDist(charge_center=(0.1, 0.1, 0.1),
                                radial_dist=[0.4, 0.4, 0.2])


@pytest.fixture
def defect_charge_info(ave_charge_density_dist):
    ave_charge_density_dist_down = AveChargeDensityDist(
        charge_center=(0.1, 0.1, 0.1), radial_dist=[0.3, 0.5, 0.2])

    return DefectChargeInfo(distance_bins=[0.0, 1.0, 2.0, 2.5],
                            band_idxs=[10],
                            charge_dists=[[ave_charge_density_dist,
                                           ave_charge_density_dist_down]],
                            ave_charge_density=0.5)


def test_charge_dist_json_roundtrip(ave_charge_density_dist, tmpdir):
    assert_json_roundtrip(ave_charge_density_dist, tmpdir)


def test_defect_charge_info_json_roundtrip(defect_charge_info, tmpdir):
    assert_json_roundtrip(defect_charge_info, tmpdir)


def test_defect_charge_info_sum_chg_dens_distribution(defect_charge_info):
    actual = defect_charge_info.sum_chg_dens_distribution(10, Spin.up)
    expected = np.array([0.4 * 1.0 ** 3,
                         0.4 * (2.0 ** 3 - 1.0 ** 3),
                         0.2 * (2.5 ** 3 - 2.0 ** 3)]) * 4 / 3 * np.pi
    np.testing.assert_array_almost_equal(actual, expected)


def test_defect_charge_info_half_charge_radius(defect_charge_info):
    actual = defect_charge_info.half_charge_radius(10, Spin.up)
    assert actual == 0.2984155182973037


def test_defect_charge_info_uniform_half_charge_radius(defect_charge_info):
    assert defect_charge_info.uniform_half_charge_radius == 0.6203504908994001


def test_defect_charge_info_repr(defect_charge_info):
    print(defect_charge_info)
    expected = """ -- defect charge info
Uniform charge radius is  0.620
Band index  Spin  Radius  Center
11          up    0.298   ( 0.100,  0.100,  0.100)
11          down  0.398   ( 0.100,  0.100,  0.100)"""
    assert defect_charge_info.__str__() == expected


def test_defect_charge_info_show_dist(defect_charge_info):
    plt = defect_charge_info.show_dist()
    plt.show()


def test_defect_charge_info_are_localized(defect_charge_info):
    actual = defect_charge_info.localized_orbitals(0.3, 1.0)
    expected = [[10], []]
    assert actual == expected

    actual = defect_charge_info.localized_orbitals(0.4, 1.0)
    expected = [[10], [10]]
    assert actual == expected

    actual = defect_charge_info.localized_orbitals(0.4, 0.398 / 0.62 - 0.001)
    expected = [[10], []]
    assert actual == expected


