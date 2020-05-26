# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List

import numpy as np
from monty.json import MSONable
from scipy.spatial import HalfspaceIntersection


@dataclass
class SingleDefectEnergy:
    name: str
    charge: int
    energy: float


@dataclass
class DefectEnergy:
    name: str
    charges: List[int]
    energies: List[float]

    def cross_points(self, ef_min, ef_max):
        large_minus_number = -1e4
        half_spaces = []
        for charge, energy in zip(self.charges, self.energies):
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
class CrossPoints:
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


