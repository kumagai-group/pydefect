# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from typing import Dict

from pymatgen import Composition


class ChemPotDiagInfo:
    def __init__(self, energies: Dict[Composition, float]):
        self.abs_energies = {c: e / c.num_atoms for c, e in energies.items()}
        self.vertex_elements = self._get_vertices()
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
            frac = [c.fractional_composition[e] for e in self.vertex_elements]
            sub = sum(f * offset for f, offset in zip(frac, self.offset_to_abs))
            result[c] = e - sub

        return result

    def _get_vertex_coords(self):

        result = [[-0.5, 0.0], [0.0, -1.0]]
        return result


