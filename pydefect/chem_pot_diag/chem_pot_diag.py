# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import string
from typing import Dict, Optional

import numpy as np
from pymatgen import Composition
from scipy.spatial.qhull import HalfspaceIntersection

from pydefect.error import PydefectError

alphabets = list(string.ascii_uppercase)


class ChemPotDiag:
    def __init__(self, energies: Dict[Composition, float]):
        self.abs_energies = {c.reduced_composition: e / c.num_atoms
                             for c, e in energies.items()}
        self.compounds = self.abs_energies.keys()
        self.vertex_elements = self._get_vertex_elements()
        self.dim = len(self.vertex_elements)
        self.offset_to_abs = self._get_offset_to_abs()
        self.rel_energies = self._get_rel_energies()
        self.vertex_coords = self._get_vertex_coords()

    def _get_vertex_elements(self):
        elements = sum(list(c.elements for c in self.compounds), [])
        return sorted(list(set(elements)))

    def _get_offset_to_abs(self):
        result = []
        for vertex_element in self.vertex_elements:
            target = Composition({vertex_element: 1.0}).reduced_composition
            candidates = filter(lambda x: x[0] == target,
                                self.abs_energies.items())
            try:
                result.append(min([x[1] for x in candidates]))
            except ValueError:
                raise NoElementEnergyError
        return result

    def _get_rel_energies(self):
        result = {}
        for c, e in self.abs_energies.items():
            sub = sum(f * offset for f, offset
                      in zip(self.frac_composition(c), self.offset_to_abs))
            result[c] = e - sub

        return result

    def frac_composition(self, c):
        return [c.fractional_composition[e] for e in self.vertex_elements]

    def _get_vertex_coords(self):
        large_minus_number = -1e5
        hs = self.get_half_space_intersection(large_minus_number)

        result = []
        for intersection in hs.intersections:
            if min(intersection) != large_minus_number:
                result.append(intersection)

        return result

    @property
    def min_rel_energies(self):
        return round(min(sum([list(l) for l in self.vertex_coords], [])) * 1.1, ndigits=4)

    def get_half_space_intersection(self, min_range):
        half_spaces = []
        for c, e in self.rel_energies.items():
            half_spaces.append(self.frac_composition(c) + [-e])
        for i in range(self.dim):
            x = [0.0] * self.dim
            x[i] = -1.0
            x.append(min_range)
            half_spaces.append(x)
        feasible_point = np.array([min_range + 1] * self.dim, dtype=float)
        hs = HalfspaceIntersection(np.array(half_spaces), feasible_point)

        return hs


class CpdPlotInfo:
    def __init__(self,
                 cpd: ChemPotDiag,
                 target: Optional[Composition] = None,
                 min_range: Optional[float] = None):
        self.cpd = cpd
        self.target = target.reduced_composition if target else None
        self.min_range = min_range or self.cpd.min_rel_energies

        self.dim = cpd.dim
        self.comp_vertices = self._get_comp_vertices(self.min_range)

    def _get_comp_vertices(self, min_range):
        hs = self.cpd.get_half_space_intersection(min_range)
        intersections = hs.intersections.tolist()
        result = {}
        for c, e in self.cpd.rel_energies.items():
            frac_comp = self.cpd.frac_composition(c)

            def on_the_composition(coord):
                diff = sum([x * y for x, y in zip(frac_comp, coord)]) - e
                return abs(diff) < 1e-8

            coords = [[round(j, ndigits=5) for j in i]
                      for i in intersections if on_the_composition(i)]

            if coords:
                result[c] = coords

        return result

    @property
    def comp_centers(self):
        return {c: np.average(np.array(v), axis=0).tolist()
                for c, v in self.comp_vertices.items()}

    def atomic_fractions(self, composition: Composition):
        return [composition.get_atomic_fraction(e)
                for e in self.cpd.vertex_elements]

    @property
    def target_vertices(self):
        label = iter(alphabets)
        vertices = self.comp_vertices[self.target]
        return {next(label): c for c in vertices if min(c) != self.min_range}


class NoElementEnergyError(PydefectError):
    pass
