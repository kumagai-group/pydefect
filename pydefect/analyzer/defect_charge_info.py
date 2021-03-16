# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from dataclasses import dataclass
from typing import List, Tuple, Sequence
import numpy as np

from monty.json import MSONable
from pymatgen import Spin
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class ChargeDist(MSONable, ToJsonFileMixIn):
    charge_center: Sequence[float]
    radial_dist: List[float]


@dataclass
class DefectChargeInfo(MSONable, ToJsonFileMixIn):
    distance_bins: List[float]  # the last value is the radius
    band_idxs: List[int]
    charge_dists: List[List[ChargeDist]]  # [band-idx, spin]
    ave_charge_density: float

    def charge_distribution(self, band_idx, spin):
        spin_idx = 0 if spin == Spin.up else 1
        band_pos = np.argwhere(np.array(self.band_idxs) == band_idx)[0][0]
        return self.charge_dists[spin_idx][band_pos].radial_dist

    def is_spin_polarized(self):
        pass

    # def is_charge_localized(self, band_idx, spin, threshold: float):
    #     _sum = 0.0
    #     for r, e in zip(self.dist(band_idx, spin), self.distance_bins):
    #         if e < self.defect_region_radius:
    #             _sum += r
    #     return _sum > threshold
    #
    # def defect_radius(self, band_idx, spin):
    #     return sum([r * e for r, e in
    #                 zip(self.dist(band_idx, spin), self.distance_bins)])

        # half_point = None
        # for i in range(hist_data[:, 1].size):
        #     if sum(hist_data[:i + 1, 1]) > 0.5:
        #         # Obtain from the calculation of the area of a trapezoid
        #         x0 = hist_data[i - 1, 0]
        #         y0 = hist_data[i, 1]
        #         y1 = hist_data[i + 1, 1]
        #         half_point = (0.5 - sum(hist_data[:i, 1]) * mesh_distance) \
        #                      * 2 / (y0 + y1) + x0
        #         break