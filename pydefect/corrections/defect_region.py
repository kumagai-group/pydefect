# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from abc import ABC, abstractmethod

import numpy as np
from monty.json import MSONable
from numpy import cross, dot
from numpy.linalg import norm
from vise.util.mix_in import ToJsonFileMixIn


def calc_max_sphere_radius(lattice_matrix) -> float:
    """Calculate Maximum radius of a sphere fitting inside the unit cell.

    Calculate three distances between two parallel planes using
    (a_i x a_j) . a_k / |a_i . a_j| """
    distances = np.zeros(3, dtype=float)
    for i in range(3):
        a_i_a_j = cross(lattice_matrix[i - 2], lattice_matrix[i - 1])
        a_k = lattice_matrix[i]
        distances[i] = abs(dot(a_i_a_j, a_k)) / norm(a_i_a_j)
    return max(distances) / 2.0


class DefectRegion(MSONable, ToJsonFileMixIn, ABC):
    @abstractmethod
    def defect_region_radius(self, lattice_matrix) -> float:
        pass


class FixedDistanceDefectRegion(DefectRegion):
    def __init__(self, radius: float):
        self.radius = radius

    def defect_region_radius(self, lattice_matrix=None) -> float:
        return self.radius


class HalfMaxFaceDistanceDefectRegion(DefectRegion):
    def __init__(self, sample_radius_ratio: float):
        self.sample_radius_ratio = sample_radius_ratio

    def defect_region_radius(self, lattice_matrix) -> float:
        sample_radius = calc_max_sphere_radius(lattice_matrix) * self.sample_radius_ratio
        return sample_radius
