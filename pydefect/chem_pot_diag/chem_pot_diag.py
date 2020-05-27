# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import string
from dataclasses import dataclass
from itertools import chain
from typing import Dict, Optional

import numpy as np
from monty.json import MSONable
from pymatgen import Composition
from scipy.spatial.qhull import HalfspaceIntersection

from pydefect.error import PydefectError
from pydefect.util.mix_in import ToJsonFileMixIn

alphabets = list(string.ascii_uppercase)


@dataclass
class ChemPotDiag(MSONable, ToJsonFileMixIn):
    energies: Dict[str, float]
    target: Optional[Composition] = None

    @property
    def abs_energies(self):
        return {Composition(c).reduced_composition: e / Composition(c).num_atoms
                for c, e in self.energies.items()}

    @property
    def compounds(self):
        return self.abs_energies.keys()

    @property
    def vertex_elements(self):
        elements = sum(list(c.elements for c in self.compounds), [])
        return sorted(list(set(elements)))

    @property
    def dim(self):
        return len(self.vertex_elements)

    @property
    def offset_to_abs(self):
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

    @property
    def rel_energies(self):
        result = {}
        for c, e in self.abs_energies.items():
            sub = sum(f * offset for f, offset
                      in zip(self.atomic_fractions(c), self.offset_to_abs))
            result[c] = e - sub

        return result

    def atomic_fractions(self, c):
        return [c.fractional_composition[e] for e in self.vertex_elements]

    @property
    def vertex_coords(self):
        large_minus_number = -1e5
        hs = self.get_half_space_intersection(large_minus_number)

        result = []
        for intersection in hs.intersections:
            if min(intersection) != large_minus_number:
                result.append(intersection)

        return result

    @property
    def min_rel_energies(self):
        return min(chain(*self.vertex_coords))

    def get_half_space_intersection(self, min_range):
        half_spaces = []
        for c, e in self.rel_energies.items():
            half_spaces.append(self.atomic_fractions(c) + [-e])
        for i in range(self.dim):
            x = [0.0] * self.dim
            x[i] = -1.0
            x.append(min_range)
            half_spaces.append(x)
        feasible_point = np.array([min_range + 1] * self.dim, dtype=float)
        hs = HalfspaceIntersection(np.array(half_spaces), feasible_point)

        return hs

    @property
    def target_vertices(self):
        label = iter(alphabets)
        fractions = self.atomic_fractions(self.target)
        energy = self.rel_energies[self.target]
        return {next(label): c for c in self.vertex_coords
                if on_composition(fractions, c, energy)}

    def abs_chem_pot_dict(self, label):
        rel_chem_pots = self.target_vertices[label]
        abs_chem_pots = [x + y for x, y in zip(rel_chem_pots, self.offset_to_abs)]
        return dict(zip(self.vertex_elements, abs_chem_pots))


class CpdPlotInfo:
    def __init__(self,
                 cpd: ChemPotDiag,
                 min_range: Optional[float] = None):
        self.cpd = cpd
        self.min_range = min_range or self.cpd.min_rel_energies * 1.1

        self.dim = cpd.dim
        self.comp_vertices = self._get_comp_vertices(self.min_range)

    def _get_comp_vertices(self, min_range):
        hs = self.cpd.get_half_space_intersection(min_range)
        intersections = hs.intersections.tolist()
        result = {}
        for c, e in self.cpd.rel_energies.items():
            def on_comp(coord):
                return on_composition(self.cpd.atomic_fractions(c), coord, e)

            coords = [[round(j, ndigits=5) for j in i]
                      for i in intersections if on_comp(i)]

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


def on_composition(atomic_fractions, coord, energy):
    diff = sum([x * y for x, y in zip(atomic_fractions, coord)]) - energy
    return abs(diff) < 1e-8


class NoElementEnergyError(PydefectError):
    pass
