# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from dataclasses import dataclass
from typing import List, Sequence

import numpy as np
from matplotlib import pyplot as plt
from monty.json import MSONable
from pymatgen import Spin
from tabulate import tabulate
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

    @property
    def uniform_half_charge_radius(self):
        return (3 * 0.5 / (4 * np.pi * self.ave_charge_density)) ** (1.0 / 3.0)

    @property
    def is_spin_polarized(self):
        return len(self.charge_dists[0]) == 2

    def half_charge_radius(self, band_idx, spin):
        charge_sum = 0.0
        spin = 0 if spin == Spin.up else 1
        idx = self.band_idxs.index(band_idx)
        dist = self.charge_dists[idx][spin].radial_dist
        for c, start, end in zip(
                dist, self.distance_bins[:-1], self.distance_bins[1:]):
            charge_sum += c
            if charge_sum > 0.5:
                return start + (end - start) * (c + 0.5 - charge_sum) / c

        raise ValueError("Radius containing 0.5 e- could not be found.")

    def __repr__(self):
        all_lines = [
            f"Uniform charge radius is {self.uniform_half_charge_radius:.3f}"]
        lines = [["band_idx", "spin", "radius"]]
        for band_idx in self.band_idxs:
            radius = self.half_charge_radius(band_idx, Spin.up)
            lines.append([band_idx, "up", round(radius, 2)])
            if self.is_spin_polarized:
                radius = self.half_charge_radius(band_idx, Spin.down)
                lines.append([band_idx, "down", round(radius, 2)])

        all_lines.append(tabulate(lines))
        return "\n".join(all_lines)

    @property
    def bins_middle_points(self):
        result = []
        for i, j in zip(self.distance_bins[:-1], self.distance_bins[1:]):
            result.append((i + j) / 2)
        return result

    def show_dist(self):
        plt.xlim([0, self.distance_bins[-1]])
        for charge_dist, _idx in zip(self.charge_dists, self.band_idxs):
            for dist, spin in zip(charge_dist, ["up", "down"]):
                plt.plot(self.bins_middle_points[:-1],
                         dist.radial_dist[:-1], label=f"{_idx} {spin}")
        ax = plt.gca()
        ax.legend(bbox_to_anchor=(0, 1.0), loc='upper left')
        return plt
