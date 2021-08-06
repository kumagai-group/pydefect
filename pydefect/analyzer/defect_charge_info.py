# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from dataclasses import dataclass
from typing import List

import numpy as np
from matplotlib import pyplot as plt
from monty.json import MSONable
from pydefect.defaults import defaults
from pydefect.util.coords import pretty_coords
from pymatgen.electronic_structure.core import Spin
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn
from vise.util.typing import Coords


@dataclass
class AveChargeDensityDist(MSONable, ToJsonFileMixIn):
    """Average charge density distribution."""
    charge_center: Coords
    radial_dist: List[float]


@dataclass
class DefectChargeInfo(MSONable, ToJsonFileMixIn):
    distance_bins: List[float]  # the last value is the radius
    band_idxs: List[int]
    charge_dists: List[List[AveChargeDensityDist]]  # [band-idx, spin]
    ave_charge_density: float

    def ave_chg_dens_distribution(self, band_idx: int, spin: Spin):
        spin_idx = 0 if spin == Spin.up else 1
        band_pos = np.argwhere(np.array(self.band_idxs) == band_idx)[0][0]
        return self.charge_dists[band_pos][spin_idx].radial_dist

    def sum_chg_dens_distribution(self, band_idx: int, spin: Spin):
        ave = self.ave_chg_dens_distribution(band_idx, spin)
        volumes = (np.array(self.distance_bins[1:]) ** 3
                   - np.array(self.distance_bins[:-1]) ** 3) * 4 / 3 * np.pi
        return list(ave * volumes)

    @property
    def uniform_half_charge_radius(self):
        return (3 * 0.5 / (4 * np.pi * self.ave_charge_density)) ** (1.0 / 3.0)

    @property
    def is_spin_polarized(self):
        return len(self.charge_dists[0]) == 2

    def half_charge_radius(self, band_idx: int, spin: Spin):
        charge_sum = 0.0
        dist = self.sum_chg_dens_distribution(band_idx, spin)
        for c, start, end in zip(
                dist, self.distance_bins[:-1], self.distance_bins[1:]):
            charge_sum += c
            if charge_sum > 0.5:
                return start + (end - start) * (c + 0.5 - charge_sum) / c

        raise ValueError("Radius containing 0.5 e- could not be found.")

    def __str__(self):
        uniform_radius = f"{self.uniform_half_charge_radius:6.3f}"
        all_lines = [" -- defect charge info",
                     f"Uniform charge radius is {uniform_radius}"]
        lines = [["Band index", "Spin", "Radius", "Center"]]
        for band_idx, c_dist in zip(self.band_idxs, self.charge_dists):
            self._add_band_info(band_idx, c_dist, lines, Spin.up)
            if self.is_spin_polarized:
                self._add_band_info(band_idx, c_dist, lines, Spin.down)

        all_lines.append(tabulate(lines, tablefmt="plain"))
        return "\n".join(all_lines)

    def _add_band_info(self, band_idx, c_dist, lines, spin):
        try:
            radius = f"{self.half_charge_radius(band_idx, spin):6.3f}"
        except ValueError:
            radius = "None"
        spin_idx = 0 if spin == Spin.up else 1
        spin_str = "up" if spin == Spin.up else "down"
        center = pretty_coords(c_dist[spin_idx].charge_center)
        lines.append([band_idx + 1, spin_str, radius, center])

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
        ax.legend(loc='upper right')
        return plt

    def localized_orbitals(self,
                           radius: float = defaults.localized_orbital_radius,
                           fraction_wrt_uniform: float
                           = defaults.localized_orbital_fraction_wrt_uniform):
        result = [[] for i in range(len(self.charge_dists[0]))]
        for charge_dist, band_idx in zip(self.charge_dists, self.band_idxs):
            for spin_idx, (dist, spin) in enumerate(zip(charge_dist,
                                                        [Spin.up, Spin.down])):
                try:
                    half_radius = self.half_charge_radius(band_idx, spin)
                    half_f = half_radius / self.uniform_half_charge_radius
                    if half_radius < radius and half_f < fraction_wrt_uniform:
                        result[spin_idx].append(band_idx)
                except ValueError:
                    pass

        return result
