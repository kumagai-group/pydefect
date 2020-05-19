# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import HalfspaceIntersection


@dataclass
class DefectEnergy:
    name: str
    charges: List[int]
    energies: List[float]
    corrections: List[float]
    shallow: List[bool]

    @property
    def corrected_energies(self):
        return [e + c for e, c in zip(self.energies, self.corrections)]

    def cross_points(self, ef_min, ef_max):
        large_minus_number = -1e4
        half_spaces = []
        for charge, energy in zip(self.charges, self.corrected_energies):
            half_spaces.append([-charge, 1, -energy])

        half_spaces.append([-1, 0, ef_min])
        half_spaces.append([1, 0, -ef_max])
        half_spaces.append([0, -1, large_minus_number])

        feasible_point = np.array([(ef_min + ef_max) / 2, -1e3])

        hs = HalfspaceIntersection(np.array(half_spaces), feasible_point)
        boundary_points = []
        inner_cross_points = []
        for intersection in hs.intersections:
            x, y = np.round(intersection, 8)
            if ef_min < x < ef_max:
                inner_cross_points.append([x, y])
            elif y > large_minus_number:
                boundary_points.append([x, y])

        all_points = boundary_points + inner_cross_points

        return (sorted(all_points, key=lambda v: v[0]),
                sorted(inner_cross_points, key=lambda v: v[0]),
                sorted(boundary_points, key=lambda v: v[0]))


@dataclass
class CrossPoints:
    inner_cross_points: List[List[float]]
    boundary_points: List[List[float]]

    @property
    def all_points(self):
        return sorted(self.boundary_points + self.inner_cross_points,
                      key=lambda v: v[0])

    @property
    def t_inner_cross_points(self):
        return np.transpose(np.array(self.inner_cross_points)).tolist()

    @property
    def t_boundary_points(self):
        return np.transpose(np.array(self.boundary_points)).tolist()


@dataclass
class DefectEnergies:
    defect_energies: List[DefectEnergy]
    vbm: float
    cbm: float
    supercell_vbm: float
    supercell_cbm: float

    def plot(self):
        for de in self.defect_energies:
            all_pts, cross_pts, boundaries = de.cross_points(self.vbm, self.cbm)
            xs, ys = np.transpose(np.array(all_pts)).tolist()
            plt.plot(xs, ys)
            xs, ys = np.transpose(np.array(cross_pts)).tolist()
            plt.scatter(xs, ys, marker="o")

        return plt

        # # Add points between x_min and x_max
        # for cp, charges in zip(tl["cross_points"], tl["charges"]):
        #     if x_min < cp[0] - self.vbm < x_max:
        #         cross_points.append([cp[0] - self.vbm, cp[1]])
        #         # need to sort the charge.
        #         charge_set.add(sorted(charges)[1])
        #         y_min, y_max = min([cp[1], y_min]), max([cp[1], y_max])
        #         ax.scatter(x=cp[0] - self.vbm,
        #                    y=cp[1],
        #                    marker='o',
        #                    color=color[i])

        # # Add the plot point at x_max
        # charge, y = \
        #     min_e_at_ef(self.defect_energies[name], x_max + self.vbm)
        # cross_points.append([x_max, y])
        # charge_set.add(charge)
        # y_min, y_max = min([y, y_min]), max([y, y_max])
        # if min(charge_set) > 0:
        #     ax.plot(x_max, y, marker="o", mec=color[i], mfc="white")
