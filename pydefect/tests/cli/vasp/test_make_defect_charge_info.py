# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import numpy as np
import pytest
from pydefect.analyzer.defect_charge_info import AveChargeDensityDist, \
    DefectChargeInfo
from pydefect.analyzer.grids import Grids
from pydefect.cli.vasp.make_defect_charge_info import \
    center_1d_periodic_quantity, make_charge_dist, make_defect_charge_info
from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp import Chgcar


def test_center_1d_periodic_quantity():
    grid_points = [0.3, 0.2, 0.1, 0.0, 0.1, 0.2, 0.3, 0.4]
    assert center_1d_periodic_quantity(grid_points) == 7
    grid_points = [0.3, 0.2, 0.1, 0.1, 0.2, 0.3, 0.4]
    assert center_1d_periodic_quantity(grid_points) == 6


@pytest.fixture
def parchg():
    struc = Structure(lattice=Lattice.cubic(10), species=["H"], coords=[[0]*3])
    data = {"total": np.array([[[0, 0], [0, 2000.0]], [[0, 0], [0, 0]]]),
            "diff": np.array([[[0, 0], [0, 1000.0]], [[0, 0], [0, 0]]])}
    # spin-up: 1.5,  spin-down: 0.5
    return Chgcar(struc, data=data)


@pytest.fixture
def charge_dist_list():
    return [AveChargeDensityDist((0.0, 0.5, 0.5), [1.5, 0.0]),
            AveChargeDensityDist((0.0, 0.5, 0.5), [0.5, 0.0])]


def test_make_charge_dist(parchg, charge_dist_list):
    grids = Grids(lattice=Lattice.cubic(10),
                  dim=(2, 2, 2),
                  distance_data=np.array(
                      [[[0., 5.], [5., 7.07106781]],
                       [[5., 7.07106781], [7.07106781, 8.66025404]]]))
    bins = np.array([0.0, 2.5, 5.0])
    actual = make_charge_dist(parchg=parchg, grids=grids, distance_bins=bins)
    # divide by 3 as there are three points in the second bin.
    assert actual == charge_dist_list


def test_make_defect_charge_info(parchg, charge_dist_list):
    grids = Grids.from_chgcar(parchg)
    actual = make_defect_charge_info([parchg], band_idxs=[1], bin_interval=2.5,
                                     grids=grids)
    expected = DefectChargeInfo([0.0, 2.5, 5.0],
                                [1],
                                [charge_dist_list],
                                0.001)
    assert actual == expected

