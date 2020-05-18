# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Dict

import numpy as np
from scipy.spatial import HalfspaceIntersection


@dataclass
class DefectEnergy:
    name: str
    num_atom_change: Dict[str, int]
    charges: List[int]
    energies: List[float]
    corrections: List[float]
    shallow: List[bool]

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

        all_points = boundary_points + inner_cross_points

        return (sorted(all_points, key=lambda v: v[0]),
                sorted(inner_cross_points, key=lambda v: v[0]))


@dataclass
class DefectEnergies:
    defect_energies: List[DefectEnergy]
    chempot: Dict[str, int]
    keywords: List[str]
    vbm: float
    cbm: float
    supercell_vbm: float
    supercell_cbm: float
