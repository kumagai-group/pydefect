# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from monty.json import MSONable
from scipy.spatial import HalfspaceIntersection


@dataclass
class DefectEnergy(MSONable):
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

        return CrossPoints(inner_cross_points, boundary_points)


@dataclass
class CrossPoints(MSONable):
    inner_cross_points: List[List[float]]
    boundary_points: List[List[float]]

    @property
    def all_sorted_points(self):
        return sorted(self.boundary_points + self.inner_cross_points,
                      key=lambda v: v[0])

    @property
    def t_all_sorted_points(self):
        return np.transpose(np.array(self.all_sorted_points)).tolist()

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
            cp = de.cross_points(self.vbm, self.cbm)
            plt.plot(*cp.t_all_sorted_points)
            plt.scatter(*cp.t_inner_cross_points, marker="o")
            plt.scatter(*cp.t_boundary_points)

        return plt
