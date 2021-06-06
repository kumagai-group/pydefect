# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import string
from dataclasses import dataclass, asdict
from itertools import product
from typing import Dict, Optional, Union, List, Set, Tuple

import numpy as np
import yaml
from monty.json import MSONable
from monty.serialization import loadfn
from pydefect.error import PydefectError
from pymatgen.core import Composition
from scipy.spatial.qhull import HalfspaceIntersection
from vise.util.logger import get_logger
from vise.util.mix_in import ToYamlFileMixIn, ToJsonFileMixIn

AtoZ = list(string.ascii_uppercase)
LargeMinusNumber = -1e5

logger = get_logger(__name__)


@dataclass
class CompositionEnergy(MSONable):
    energy: float
    source: str = None


class CompositionEnergies(ToYamlFileMixIn, dict):
    """
    keys: Composition
    values: CompositionEnergy
    """
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

    @classmethod
    def from_dict(cls, d: Dict[str, float]):
        ce = {Composition(k): CompositionEnergy(v) for k, v in d.items()}
        return cls(ce)

    @property
    def elements(self):
        result = set()
        for c in self:
            result.update(set([str(e) for e in c.elements]))
        return sorted(result)

    @property
    def std_rel_energies(self) -> Tuple["StandardEnergies", "RelativeEnergies"]:
        std, rel = StandardEnergies(), RelativeEnergies()
        abs_energies_per_atom = {k.reduced_formula: v.energy / k.num_atoms
                                 for k, v in self.items()}
        std_energies_list = []
        for vertex_element in self.elements:
            # This target is needed as some reduced formulas shows molecule
            # ones such as H2 and O2.
            reduced_formula = Composition({vertex_element: 1.0}).reduced_formula
            candidates = filter(lambda x: x[0] == reduced_formula,
                                abs_energies_per_atom.items())
            try:
                min_abs_energy = min([abs_energy_per_atom[1]
                                      for abs_energy_per_atom in candidates])
            except ValueError:
                print(f"Element {vertex_element} does not exist in "
                      f"CompositionEnergies.")
                raise NoElementEnergyError
            std[vertex_element] = min_abs_energy
            std_energies_list.append(min_abs_energy)

        for formula, abs_energy_per_atom in abs_energies_per_atom.items():
            if Composition(formula).is_element:
                continue
            frac = atomic_fractions(formula, self.elements)
            offset = sum([a * b for a, b in zip(frac, std_energies_list)])
            rel[formula] = abs_energy_per_atom - offset

        return std, rel


class CpdAbstractEnergies(ToYamlFileMixIn, dict):
    """
    keys: str (composition name)
    values: float (energy per atom)
    """
    def to_yaml(self) -> str:
        return yaml.dump(dict(self))

    @classmethod
    def from_yaml(cls, filename: str = None):
        return cls(loadfn(filename or cls._yaml_filename()))


class StandardEnergies(CpdAbstractEnergies):
    pass


def atomic_fractions(comp: Union[Composition, str], elements: List[str]
                     ) -> List[float]:
    return [Composition(comp).fractional_composition[e] for e in elements]


def comp_to_element_set(comp: Union[Composition, str]) -> Set[str]:
    return {str(e) for e in Composition(comp).elements}


def target_element_chem_pot(comp: Union[Composition, str],
                            energy_per_atom: float,
                            target_element: str,
                            other_elem_chem_pot: Dict[str, float]) -> float:
    assert comp_to_element_set(comp) \
           <= set(other_elem_chem_pot) | {target_element}
    other_element_val = 0.0
    for element, frac in Composition(comp).fractional_composition.items():
        if target_element == str(element):
            target_frac = frac
            continue
        other_element_val += other_elem_chem_pot[str(element)] * frac
    return (energy_per_atom - other_element_val) / target_frac


class RelativeEnergies(CpdAbstractEnergies):
    #TODO: add stability
    @property
    def all_element_set(self) -> Set[str]:
        return set().union(*[comp_to_element_set(c) for c in self])

    def host_composition_energies(self, elements: List[str]
                                  ) -> Dict[str, float]:
        return {formula: energy for formula, energy in self.items()
                if comp_to_element_set(formula).issubset(elements)}

    def comp_energies_with_element(self,
                                   element: str) -> Dict[str, float]:
        return {formula: energy for formula, energy in self.items()
                if element in comp_to_element_set(formula)}

    def impurity_chem_pot(self, impurity_element: str,
                          host_elements_chem_pot: Dict[str, float]
                          ) -> Tuple[float, str]:
        impurity_chem_pot = {impurity_element: 0.0}
        comp_energies = self.comp_energies_with_element(impurity_element)

        for formula, energy_per_atom in comp_energies.items():
            impurity_chem_pot[formula] = target_element_chem_pot(
                formula, energy_per_atom, impurity_element,
                host_elements_chem_pot)

        competing_phase_formula = min(impurity_chem_pot,
                                      key=impurity_chem_pot.get)
        chem_pot = impurity_chem_pot[competing_phase_formula]
        return chem_pot, competing_phase_formula


class ChemPotDiagMaker:
    def __init__(self,
                 relative_energies: RelativeEnergies,
                 elements: List[str],
                 target: str = None):
        self.relative_energies = relative_energies
        self.host_composition_energies = \
            relative_energies.host_composition_energies(elements)
        self.elements = elements
        self.impurity_elements = \
            relative_energies.all_element_set.difference(elements)
        self.dim = len(elements)

        if target:
            try:
                assert target in relative_energies.keys()
            except AssertionError:
                print(f"Target {target} is not in relative energy compounds.")
                raise
        self.target = target

    def _calc_vertices(self):
        half_spaces = []

        for formula, energy in self.host_composition_energies.items():
            half_spaces.append(
                atomic_fractions(formula, self.elements) + [-energy])

        for i in range(self.dim):
            upper_boundary, lower_boundary = [0.0] * self.dim, [0.0] * self.dim
            upper_boundary[i], lower_boundary[i] = 1.0, -1.0

            upper_boundary.append(0.0)
            lower_boundary.append(LargeMinusNumber)
            half_spaces.extend([upper_boundary, lower_boundary])

        feasible_point = np.array([LargeMinusNumber + 1.0] * self.dim,
                                  dtype=float)
        hs = HalfspaceIntersection(np.array(half_spaces), feasible_point)
        self.vertices: List[List[float]] = hs.intersections.tolist()

    def _is_composition_involved(self,
                                 coord: List[float],
                                 composition: str,
                                 energy: float) -> bool:
        atom_frac = atomic_fractions(composition, self.elements)
        diff = sum([x * y for x, y in zip(atom_frac, coord)]) - energy
        return abs(diff) < 1e-8

    def _min_energy_range(self, mul: float = 1.1):
        vertex_values = [x for x in sum(self.vertices, [])
                         if x != LargeMinusNumber]
        return min(vertex_values) * mul

    @property
    def chem_pot_diag(self):
        self._calc_vertices()
        elem_energies = {element: 0.0 for element in self.elements}
        host_energies = dict(**self.host_composition_energies, **elem_energies)

        polygons = {}
        for comp, energy in host_energies.items():
            vertices = []
            for coord in self.vertices:
                if self._is_composition_involved(coord, comp, energy):
                    vertex = [round(c, ndigits=5) if c != LargeMinusNumber
                              else self._min_energy_range() for c in coord]
                    vertices.append(vertex)
            if vertices:
                polygons[comp] = vertices

        vertices = None
        if self.target:
            target_vertices = []
            for coord in polygons[self.target]:
                competing_phases = []
                for comp, vertices in polygons.items():
                    if comp == self.target:
                        continue
                    if coord in vertices:
                        competing_phases.append(comp)

                impurity_phases = []
                host_chem_pots = dict(zip(self.elements, coord))
                impurity_chem_pots = {}
                for i_element in self.impurity_elements:
                    i_chem_pot, i_phase = \
                        self.relative_energies.impurity_chem_pot(
                            i_element, host_chem_pots)
                    impurity_chem_pots[i_element] = i_chem_pot
                    impurity_phases.append(i_phase)
                chem_pot = dict(**host_chem_pots, **impurity_chem_pots)
                target_vertices.append(TargetVertex(chem_pot,
                                                    competing_phases,
                                                    impurity_phases))

            AtoZZ = product([""] + AtoZ, AtoZ)
            vertices = {"".join(next(AtoZZ)): v for v in target_vertices}

        return ChemPotDiag(vertex_elements=self.elements,
                           polygons=polygons,
                           target=self.target,
                           target_vertices_dict=vertices)


@dataclass
class TargetVertex(MSONable):
    chem_pot: Dict[str, float]
    competing_phases: Optional[List[str]] = None
    impurity_phases: Optional[List[str]] = None


@dataclass
class TargetVertices(ToYamlFileMixIn):
    target: str
    vertices: Dict[str, TargetVertex]

    @property
    def chem_pots(self) -> Dict[str, Dict[str, float]]:
        return {k: v.chem_pot for k, v in self.vertices.items()}

    def to_yaml(self) -> str:
        header = f"target: {self.target}"
        d = {k: asdict(v) for k, v in self.vertices.items()}
        return "\n".join([header, yaml.dump(d)])

    @classmethod
    def from_yaml(cls, filename: str = None):
        name = filename or cls._yaml_filename()
        d = loadfn(name)
        target = d.pop("target")
        vertices = {k: TargetVertex(**v) for k, v in d.items()}
        return cls(target=target, vertices=vertices)


@dataclass
class ChemPotDiag(MSONable, ToJsonFileMixIn):
    vertex_elements: List[str]
    polygons: Dict[str, List[List[float]]]
    target: str = None
    target_vertices_dict: Dict[str, TargetVertex] = None

    @property
    def target_coords(self):
        if self.target_vertices_dict:
            return {k: [v.chem_pot[e] for e in self.vertex_elements]
                    for k, v in self.target_vertices_dict.items()}

    @property
    def min_range(self):
        return np.min(sum(self.polygons.values(), []))

    @property
    def dim(self):
        return len(self.vertex_elements)

    @property
    def comp_centers(self):
        _polygons = {}
        return {c: np.average(np.array(v), axis=0).tolist()
                for c, v in self.polygons.items()}

    def atomic_fractions(self, composition: str):
        return [Composition(composition).get_atomic_fraction(e)
                for e in self.vertex_elements]

    @property
    def chemical_system(self) -> str:
        return "-".join([el for el in self.vertex_elements])

    @property
    def to_target_vertices(self):
        if self.target and self.target_vertices_dict:
            return TargetVertices(self.target, self.target_vertices_dict)
        logger.warning("Need to set target and target_vertices.")
        raise ValueError


class NoElementEnergyError(PydefectError):
    pass
