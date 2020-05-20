# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import Dict

import numpy as np
from pymatgen import Composition
from scipy.spatial.qhull import HalfspaceIntersection


class ChemPotDiag:
    def __init__(self, energies: Dict[Composition, float]):
        self.abs_energies = {c: e / c.num_atoms for c, e in energies.items()}
        self.vertex_elements = self._get_vertices()
        self.dim = len(self.vertex_elements)
        self.offset_to_abs = self._get_offset_to_abs()
        self.rel_energies = self._get_rel_energies()
        self.vertex_coords = self._get_vertex_coords()

    def _get_vertices(self):
        elements = sum(list(c.elements for c in self.abs_energies.keys()), [])
        return sorted(list(set(elements)))

    def _get_offset_to_abs(self):
        result = []
        for vertex_element in self.vertex_elements:
            target = Composition(str(vertex_element)).reduced_formula
            candidates = filter(lambda x: x[0].reduced_formula == target,
                                self.abs_energies.items())
            result.append(min([x[1] for x in candidates]))
        return result

    def _get_rel_energies(self):
        result = {}
        for c, e in self.abs_energies.items():
            sub = sum(f * offset for f, offset
                      in zip(self._frac_composition(c), self.offset_to_abs))
            result[c] = e - sub

        return result

    def _frac_composition(self, c):
        return [c.fractional_composition[e] for e in self.vertex_elements]

    def _get_vertex_coords(self):
        large_minus_number = -1e5
        half_spaces = []
        for c, e in self.rel_energies.items():
            half_spaces.append(self._frac_composition(c) + [-e])

        for i in range(self.dim):
            x = [0.0] * self.dim
            x[i] = -1.0
            x.append(large_minus_number)
            half_spaces.append(x)

        feasible_point = np.array([large_minus_number + 1] * self.dim)
        hs = HalfspaceIntersection(np.array(half_spaces), feasible_point)

        result = []
        for intersection in hs.intersections:
            if min(intersection) != large_minus_number:
                result.append(intersection)

        return result


