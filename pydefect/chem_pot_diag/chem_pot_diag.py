# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import string
from copy import deepcopy
from dataclasses import dataclass, InitVar
from itertools import chain
from pathlib import Path
from typing import Dict, Optional, Union, List, Set

import numpy as np
import yaml
from monty.serialization import loadfn
from pydefect.error import PydefectError
from pydefect.util.error_classes import CpdNotSupportedError
from pymatgen import Composition, Element
from scipy.spatial.qhull import HalfspaceIntersection

alphabets = list(string.ascii_uppercase)


@dataclass(frozen=True)
class CompositionEnergy:
    composition: Composition
    energy: float
    source: str

    @property
    def abs_energy_per_atom(self):
        return self.energy / self.composition.num_atoms


@dataclass
class ChemPotDiag:
    comp_energies: Set[CompositionEnergy]
    target: InitVar[Union[Composition, dict]]
    vertex_elements: InitVar[Optional[List[Element]]] = None

    def __post_init__(self, target, vertex_elements):
        self.target: Composition = target if isinstance(target, Composition) \
            else Composition.from_dict(target)
        self.vertex_elements = vertex_elements or sorted(self.target.elements)

    def to_yaml(self, filename: str = "cpd.yaml") -> None:
        d = {"target": str(self.target).replace(" ", ""),
             "vertex_elements": [str(e) for e in self.vertex_elements]}
        for ce in self.comp_energies:
            key = str(ce.composition.iupac_formula).replace(" ", "")
            val = {"energy": ce.energy, "source": ce.source}
            d[key] = val
        Path(filename).write_text(yaml.dump(d))

    @classmethod
    def from_yaml(cls, filename: str = "cpd.yaml") -> "ChemPotDiag":
        d = loadfn(filename)
        target = d.pop("target")
        vertex_elements = [Element(e) for e in d.pop("vertex_elements")]
        composition_energies = set()
        for k, v in d.items():
            composition_energies.add(
                CompositionEnergy(Composition(k), v["energy"], v["source"]))
        return cls(composition_energies, target, vertex_elements)

    @property
    def abs_energies_per_atom(self) -> Dict[Composition, float]:
        return {ce.composition.reduced_composition: ce.abs_energy_per_atom
                for ce in self.comp_energies}

    @property
    def all_compounds(self) -> List[Composition]:
        return list(self.abs_energies_per_atom.keys())

    @property
    def impurity_elements(self):
        result = set()
        for c in self.all_compounds:
            result.update(set(c.elements))
        result.difference_update(set(self.target.elements))
        return list(result)

    @property
    def dim(self):
        return len(self.vertex_elements)

    @property
    def host_ele_abs_energies_per_atom(self):
        result = {}
        for k, v in self.abs_energies_per_atom.items():
            if set(k.elements).issubset(set(self.vertex_elements)):
                result[k] = v
        return result

    @property
    def offset_to_abs(self):
        result = []
        for vertex_element in self.vertex_elements:
            target = Composition({vertex_element: 1.0}).reduced_composition
            candidates = filter(lambda x: x[0] == target,
                                self.host_ele_abs_energies_per_atom.items())
            try:
                result.append(min([x[1] for x in candidates]))
            except ValueError:
                raise NoElementEnergyError
        return result

    @property
    def rel_energies(self):
        result = {}
        for c, e in self.host_ele_abs_energies_per_atom.items():
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
    def lowest_relative_energy(self):
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

    def host_abs_chem_pot_dict(self, label) -> Dict[Element, float]:
        rel_chem_pots = self.target_vertices[label]
        abs_chem_pots = [x + y for x, y in zip(rel_chem_pots, self.offset_to_abs)]
        return dict(zip(self.vertex_elements, abs_chem_pots))

    def abs_chem_pot_dict(self, label) -> Dict[Element, float]:
        result = deepcopy(self.host_abs_chem_pot_dict(label))
        for e in self.impurity_elements:
            _, result[e] = self.impurity_abs_energy(e, label)
        return result

    def impurity_abs_energy(self, element: Element, label: str):
        comp_set = set(self.vertex_elements) | {element}
        competing_comp_e = None
        y = float("inf")
        for ce in self.comp_energies:
            if element in set(ce.composition.elements):
                if set(ce.composition.elements).issubset(comp_set) is False:
                    raise CpdNotSupportedError(
                        "Other element(s) than host elements exists.")
                abs_chem_pot = self.host_abs_chem_pot_dict(label)
                mu = ce.energy
                comp_d = ce.composition.as_dict()
                for ve in self.vertex_elements:
                    mu -= comp_d[str(ve)] * abs_chem_pot[ve]
                mu /= comp_d[str(element)]
                if mu < y:
                    competing_comp_e, y = ce, mu
        if competing_comp_e is None:
            raise CpdNotSupportedError(f"No compounds with element {element}.")
        return competing_comp_e, y


class CpdPlotInfo:
    def __init__(self,
                 cpd: ChemPotDiag,
                 min_range: Optional[float] = None):
        self.cpd = cpd
        self.min_range = min_range or self.cpd.lowest_relative_energy * 1.1

        self.dim = cpd.dim
        self.comp_vertices = self._get_comp_vertices(self.min_range)

    def _get_comp_vertices(self, min_range
                           ) -> Dict[Composition, List[List[float]]]:
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


def on_composition(atomic_fractions, coord, energy) -> bool:
    diff = sum([x * y for x, y in zip(atomic_fractions, coord)]) - energy
    return abs(diff) < 1e-8


def replace_comp_energy(chem_pot_diag: ChemPotDiag,
                        replaced_comp_energies: Set[CompositionEnergy]):
    new_comp_energies = set()
    for ce in chem_pot_diag.comp_energies:
        for replaced_comp_energy in replaced_comp_energies:
            if (ce.composition.reduced_composition
                    == replaced_comp_energy.composition.reduced_composition):
                new_comp_energies.add(replaced_comp_energy)
                break
        else:
            new_comp_energies.add(ce)
    chem_pot_diag.comp_energies = new_comp_energies


class NoElementEnergyError(PydefectError):
    pass
