# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import string
from copy import deepcopy
from dataclasses import dataclass, InitVar
from itertools import chain
from pathlib import Path
from typing import Dict, Optional, Union, List

import numpy as np
import pandas as pd
import yaml
from monty.json import MSONable
from monty.serialization import loadfn
from pydefect.error import PydefectError
from pydefect.util.error_classes import CpdNotSupportedError
from pymatgen import Composition, Element
from scipy.spatial.qhull import HalfspaceIntersection, QhullError
from tabulate import tabulate

alphabets = list(string.ascii_uppercase) + \
            ["A"+i for i in string.ascii_uppercase] + \
            ["B"+i for i in string.ascii_uppercase] + \
            ["C"+i for i in string.ascii_uppercase]


@dataclass
class CompositionEnergy(MSONable):
    composition: InitVar[Union[Composition, dict]]
    energy: float
    source: str

    def __post_init__(self, composition):
        self.composition: Composition = composition if isinstance(composition, Composition) \
            else Composition.from_dict(composition)

    @property
    def abs_energy_per_atom(self):
        return self.energy / self.composition.num_atoms


@dataclass
class ChemPotDiag(MSONable):
    comp_energies: List[CompositionEnergy]
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
        composition_energies = []
        for k, v in d.items():
            composition_energies.append(
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
    def vertex_elements_abs_energies_per_atom(self):
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
                                self.vertex_elements_abs_energies_per_atom.items())
            try:
                result.append(min([x[1] for x in candidates]))
            except ValueError:
                raise NoElementEnergyError
        return result

    @property
    def rel_energies(self):
        result = {}
        for c, e in self.vertex_elements_abs_energies_per_atom.items():
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

    def get_half_space_intersection(self, min_range, feasible_dist=1.0):
        half_spaces = []
        for c, e in self.rel_energies.items():
            half_spaces.append(self.atomic_fractions(c) + [-e])
        for i in range(self.dim):
            x = [0.0] * self.dim
            x[i] = -1.0
            x.append(min_range)
            half_spaces.append(x)
        # When used in vertex_coords  min_range is large_minum_number
        feasible_point = np.array([min_range + feasible_dist] * self.dim, dtype=float)
        hs = HalfspaceIntersection(np.array(half_spaces), feasible_point)

        return hs

    @property
    def target_vertices(self):
        label = iter(alphabets)
        fractions = self.atomic_fractions(self.target)
        energy = self.rel_energies[self.target]
        return {next(label): c for c in self.vertex_coords
                if on_composition(fractions, c, energy)}

    @property
    def target_vertex_list_dataframe(self):
        index = []
        result = []
        for k, v in self.target_vertices.items():
            index.append(k)
            result.append([])
            for k2, v2 in zip(self.vertex_elements, v):
                result[-1].append(round(v2, 3))

            for ie in self.impurity_elements:
                competing_comp_for_impurity, _ = self.impurity_abs_energy(ie, k)
                comp_name = str(competing_comp_for_impurity.composition)
                result[-1].append(comp_name)

        columns = [f"mu_{e}" for e in self.vertex_elements]
        columns += [f"Phase for {e}" for e in self.impurity_elements]
        return pd.DataFrame(result, index=index, columns=columns)

    def __str__(self):
        return tabulate(self.target_vertex_list_dataframe, headers='keys', tablefmt='psql')

    def target_abs_chem_pot_dict(self, label) -> Dict[Element, float]:
        rel_chem_pots = self.target_vertices[label]
        abs_chem_pots = [x + y for x, y in zip(rel_chem_pots, self.offset_to_abs)]
        return dict(zip(self.vertex_elements, abs_chem_pots))

    def abs_chem_pot_dict(self, label) -> Dict[Element, float]:
        result = deepcopy(self.target_abs_chem_pot_dict(label))
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
                abs_chem_pot = self.target_abs_chem_pot_dict(label)
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

    @property
    def target_formation_energy(self):
        return self.rel_energies[self.target]


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
        try:
            hs = self.cpd.get_half_space_intersection(min_range)
        except QhullError:
            hs = self.cpd.get_half_space_intersection(min_range,
                                                      feasible_dist=0.01)

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
                        replaced_comp_energies: List[CompositionEnergy]):
    new_comp_energies = []
    for ce in chem_pot_diag.comp_energies:
        for replaced_comp_energy in replaced_comp_energies:
            if (ce.composition.reduced_composition
                    == replaced_comp_energy.composition.reduced_composition):
                new_comp_energies.append(replaced_comp_energy)
                break
        else:
            new_comp_energies.append(ce)
    chem_pot_diag.comp_energies = new_comp_energies


class NoElementEnergyError(PydefectError):
    pass
