# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import string
from copy import copy
from dataclasses import dataclass
from itertools import product
from typing import Dict, Optional, Union, List, Set

import numpy as np
import yaml
from monty.json import MSONable
from monty.serialization import loadfn
from pydefect.error import PydefectError
from pymatgen.core import Composition
from scipy.spatial.qhull import HalfspaceIntersection
from vise.util.mix_in import ToYamlFileMixIn

AtoZ = list(string.ascii_uppercase)
LargeMinusNumber = -1e5


@dataclass
class CompositionEnergy(MSONable):
    energy: float
    source: str = None


class CompositionEnergies(ToYamlFileMixIn, dict):
    def to_yaml(self) -> str:
        d = {}
        for k, v in self.items():
            key = str(k.iupac_formula).replace(" ", "")
            val = {"energy": v.energy, "source": str(v.source)}
            d[key] = val
        return yaml.dump(d)

    @classmethod
    def from_yaml(cls, filename: str = None):
        name = filename or cls._yaml_filename()
        d = loadfn(name)
        composition_energies = {}
        for k, v in d.items():
            source = v.get('source', None)
            key = Composition(k)
            composition_energies[key] = CompositionEnergy(v["energy"], source)
        return cls(composition_energies)

    @property
    def elements(self):
        result = set()
        for c in self:
            result.update(set([str(e) for e in c.elements]))
        return sorted(result)

    @property
    def std_rel_energies(self):
        std, rel = StandardEnergies(), RelativeEnergies()
        abs_energies_per_atom = {k.reduced_formula: v.energy / k.num_atoms
                                 for k, v in self.items()}
        ref_energy_list = []
        for vertex_element in self.elements:
            # This target is needed as some reduced formulas shows molecule
            # ones such as H2 and O2.
            target = Composition({vertex_element: 1.0}).reduced_formula
            candidates = filter(lambda x: x[0] == target,
                                abs_energies_per_atom.items())
            try:
                min_abs_energy = min([x[1] for x in candidates])
            except ValueError:
                raise NoElementEnergyError
            std[vertex_element] = min_abs_energy
            ref_energy_list.append(min_abs_energy)

        for formula, e in abs_energies_per_atom.items():
            comp = Composition(formula)
            if comp.is_element:
                continue
            frac = [comp.fractional_composition[e] for e in self.elements]
            offset = sum([a * b for a, b in zip(frac, ref_energy_list)])
            rel[formula] = e - offset

        return std, rel


class CpdAbstractEnergies(ToYamlFileMixIn, dict):
    def to_yaml(self) -> str:
        return yaml.dump(dict(self))

    @classmethod
    def from_yaml(cls, filename: str = None):
        return cls(loadfn(filename or cls._yaml_filename()))


class StandardEnergies(CpdAbstractEnergies):
    pass


def atomic_fractions(comp: Union[Composition, str], elements: List[str]):
    return [Composition(comp).fractional_composition[e] for e in elements]


def comp_to_element_set(comp: str) -> Set[str]:
    return {str(e) for e in Composition(comp).elements}


def calc_chem_pot(comp, energy_per_atom, target_elem, other_elem_chem_pot):
    x = 0.0
    fracs = Composition(comp).fractional_composition
    for element, frac in fracs.items():
        if target_elem == str(element):
            continue
        x += other_elem_chem_pot[str(element)] * frac
    return (energy_per_atom - x) / fracs[target_elem]


class RelativeEnergies(CpdAbstractEnergies):
    #TODO: add stability
    @property
    def all_element_set(self) -> Set[str]:
        return set().union(*[comp_to_element_set(c) for c in self])

    def related_comp_energies(self, elements: List[str]) -> Dict[str, float]:
        return {c: e for c, e in self.items()
                if comp_to_element_set(c).issubset(elements)}

    def comp_energies_w_element(self, element: str) -> Dict[str, float]:
        return {c: e for c, e in self.items()
                if element in comp_to_element_set(c)}

    def impurity_info(self,
                      impurity: str,
                      host_chem_pot: Dict[str, float]) -> "ImpurityInfo":
        candidates = {impurity: 0.0}
        for c, e in self.comp_energies_w_element(impurity).items():
            candidates[c] = calc_chem_pot(c, e, impurity, host_chem_pot)

        min_key = min(candidates, key=candidates.get)
        return ImpurityInfo(candidates[min_key], min_key)


class ChemPotDiagMaker:
    def __init__(self,
                 rel_energies: RelativeEnergies,
                 elements: List[str],
                 target: str = None):
        self.rel_energies = rel_energies
        self.comp_energies = rel_energies.related_comp_energies(elements)
        self.elements = elements
        self.impurities = rel_energies.all_element_set.difference(elements)
        self.dim = len(elements)
        if target:
            try:
                assert target in rel_energies.keys()
            except AssertionError:
                print(f"Target {target} is not in relative energy compounds.")
        self.target = target

    def _set_vertices(self):
        half_spaces = []

        for c, e in self.comp_energies.items():
            half_spaces.append(atomic_fractions(c, self.elements) + [-e])

        for i in range(self.dim):
            upper_boundary, lower_boundary = [0.0] * self.dim, [0.0] * self.dim
            upper_boundary[i], lower_boundary[i] = 1.0, -1.0

            upper_boundary.append(0.0)
            lower_boundary.append(LargeMinusNumber)
            half_spaces.extend([upper_boundary, lower_boundary])

        # When used in vertex_coords min_range is large_minus_number
        feasible_point = np.array([LargeMinusNumber + 1.0] * self.dim,
                                  dtype=float)
        hs = HalfspaceIntersection(np.array(half_spaces), feasible_point)
        self.vertices = hs.intersections.tolist()

    def _on_comp(self, coord, composition, energy):
        atom_frac = atomic_fractions(composition, self.elements)
        return on_composition(atom_frac, coord, energy)

    def _min_value(self, mul: float = 1.1):
        candidates = [x for x in sum(self.vertices, []) if x != LargeMinusNumber]
        return min(candidates) * mul

    def _polygons(self):
        self._set_vertices()
        result = {}
        comp_energies = copy(self.comp_energies)

        if self.target:
            target_coords = [[round(j, ndigits=5) for j in i]
                              for i in self.vertices
                              if self._on_comp(i, self.target,
                                               self.comp_energies[self.target])]
            competing_phases = [[] for i in target_coords]

        min_val = self._min_value()

        # elements
        for _idx, element in enumerate(self.elements):
            coords = [[round(j, ndigits=5) for j in i]
                      for i in self.vertices if round(i[_idx], ndigits=5) == 0.0]
            coords = [[c if c != LargeMinusNumber else min_val for c in cc] for cc in coords]
            result[element] = coords

            if self.target:
                for c in coords:
                    if c in target_coords:
                        _idx = target_coords.index(c)
                        competing_phases[_idx].append(element)

        # compounds
        for comp, energy in comp_energies.items():
            coords = [[round(j, ndigits=5) for j in i]
                      for i in self.vertices if self._on_comp(i, comp, energy)]
            if coords:
                result[comp] = coords
                if self.target and comp != self.target:
                    for c in coords:
                        if c in target_coords:
                            _idx = target_coords.index(c)
                            competing_phases[_idx].append(comp)

        x = []
        if self.target:
            for v, w in zip(target_coords, competing_phases):
                host_chem_pot = dict(zip(self.elements, v))
                i = {imp: self.rel_energies.impurity_info(imp, host_chem_pot)
                     for imp in self.impurities}
                x.append(TargetVertex(v, w, i))
        return result, x or None

    @property
    def chem_pot_diag(self):
        AtoZZ = product([""] + AtoZ, AtoZ)
        polygons, target_vertices = self._polygons()
        if self.target:
            target_vertices = TargetVertices(self.target, {"".join(next(AtoZZ)): v for v in target_vertices})
        else:
            target_vertices = None
        return ChemPotDiag(vertex_elements=self.elements,
                           polygons=polygons,
                           target_vertices=target_vertices)


@dataclass
class ImpurityInfo:
    rel_chem_pot: float
    competing_phase: str


@dataclass
class TargetVertex(MSONable):
    coords: List[float]
    competing_phases: Optional[List[str]]
    impurity_info: Optional[Dict[str, ImpurityInfo]]

    @property
    def impurity_chem_pot(self):
        return {k: v.rel_chem_pot for k, v in self.impurity_info.items()}


@dataclass
class TargetVertices:
    target: str
    vertices: Dict[str, TargetVertex]


@dataclass
class ChemPotDiag(MSONable):
    vertex_elements: List[str]
    polygons: Dict[str, List[List[float]]]
    target_vertices: TargetVertices = None

    @property
    def min_value(self):
        return np.min(sum(self.polygons.values(), []))

    @property
    def target(self):
        if self.target_vertices:
            return self.target_vertices.target
        return

    @property
    def dim(self):
        return len(self.vertex_elements)

    @property
    def comp_centers(self):
        _polygons = {}
        return {c: np.average(np.array(v), axis=0).tolist()
                for c, v in self.polygons.items()}

    def atomic_fractions(self, composition: Composition):
        return [composition.get_atomic_fraction(e)
                for e in self.vertex_elements]


def on_composition(atomic_fracs, coord, energy) -> bool:
    diff = sum([x * y for x, y in zip(atomic_fracs, coord)]) - energy
    return abs(diff) < 1e-8


# def replace_comp_energy(chem_pot_diag: ChemPotDiag,
#                         replaced_comp_energies: List[CompositionEnergy]):
#     new_comp_energies = []
#     for ce in chem_pot_diag.comp_energies:
#         for replaced_comp_energy in replaced_comp_energies:
#             if (ce.composition.reduced_composition
#                     == replaced_comp_energy.composition.reduced_composition):
#                 new_comp_energies.append(replaced_comp_energy)
#                 break
#         else:
#             new_comp_energies.append(ce)
#     chem_pot_diag.comp_energies = new_comp_energies


class NoElementEnergyError(PydefectError):
    pass
