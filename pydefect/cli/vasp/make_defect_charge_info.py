# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from itertools import product
from typing import List
import numpy as np
from pydefect.analyzer.defect_charge_info import DefectChargeInfo, ChargeDist
from pydefect.analyzer.grids import Grids
from pydefect.cli.vasp.make_efnv_correction import calc_max_sphere_radius
from pymatgen import Spin
from pymatgen.io.vasp import Chgcar


def center_1d_periodic_quantity(grid_points: List[float]):
    num_grid_pts = len(grid_points)
    moments = []
    for i in range(num_grid_pts):
        moment = 0.0
        if num_grid_pts % 2 == 0:
            for j in range(1, int(num_grid_pts / 2 + 1)):
                moment += grid_points[i - j] * j
            for j in range(1, int(num_grid_pts / 2)):
                k = (i + j - num_grid_pts) % num_grid_pts
                moment += grid_points[k] * j
        else:
            a = int(np.floor(num_grid_pts / 2))
            for j in range(1, a + 1):
                moment += grid_points[i - j] * j
            for j in range(1, a + 1):
                k = (i + j - num_grid_pts) % num_grid_pts
                moment += grid_points[k] * j

        moments.append(moment)

    return np.nanargmin(moments)


def make_spin_charges(parchg: Chgcar) -> List[Chgcar]:
    result = [Chgcar(parchg.structure, {"total": parchg.spin_data[Spin.up]})]
    if "diff" in parchg.data:
        result.append(
            Chgcar(parchg.structure, {"total": parchg.spin_data[Spin.down]}))
    return result


def make_charge_dist(parchg: Chgcar, grids: Grids, distance_bins: np.ndarray):
    assert parchg.structure.lattice == grids.lattice
    assert parchg.dim == grids.dim

    spin_charges = make_spin_charges(parchg)

    dists = []
    defect_center_idxs = []
    for chg in spin_charges:
        center = [center_1d_periodic_quantity(
            chg.get_average_along_axis(i)) for i in [0, 1, 2]]
        defect_center_idxs.append(np.array(center))
        data = chg.data["total"] / grids.lattice.volume
        dists.append(grids.spherical_dist(data, center, distance_bins))

    result = []
    for c, d in zip(defect_center_idxs, dists):
        result.append(ChargeDist((c / grids.dim).tolist(), d))

    return result


def make_defect_charge_info(parchgs: List[Chgcar],
                            band_idxs: List[int],
                            bin_interval: float):
    grids = Grids.from_chgcar(parchgs[0])
    radius = calc_max_sphere_radius(parchgs[0].structure.lattice.matrix)
    num_bins = int(np.ceil(radius / bin_interval))
    distance_bins = np.array([bin_interval * i for i in range(num_bins)] + [radius])
    ave_charge_density = 1.0 / parchgs[0].structure.volume
    charge_dists = []
    for parchg in parchgs:
        charge_dists.append(make_charge_dist(parchg, grids, distance_bins))
    return DefectChargeInfo(distance_bins.tolist(), band_idxs, charge_dists, ave_charge_density)

