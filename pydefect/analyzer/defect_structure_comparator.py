# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from dataclasses import dataclass
from typing import List, Tuple, Dict

import numpy as np
from monty.json import MSONable
from pydefect.util.structure_tools import Distances
from pymatgen import IStructure


class DefectStructureComparator:
    def __init__(self,
                 defect_structure: IStructure,
                 perfect_structure: IStructure,
                 dist_tol: float = None):
        """
        Atoms in the final structure are shifted such that the farthest atom
        from the defect is placed at the same place with that in the perfect
        supercell.
        """
        self._defect_structure = defect_structure
        self._perfect_structure = perfect_structure
        self.dist_tol = dist_tol
        self.p_to_d = self.make_p_to_d()
        self.d_to_p = self.make_d_to_p()

    @property
    def atom_mapping(self):
        return {d: p for d, p in enumerate(self.d_to_p)
                if d not in self.inserted_indices}

    def make_p_to_d(self):
        result = []
        for site in self._perfect_structure:
            distances = Distances(self._defect_structure,
                                  site.frac_coords,
                                  self.dist_tol)
            result.append(distances.atom_idx_at_center(specie=site.specie))
        return result

    def make_d_to_p(self):
        result = []
        for site in self._defect_structure:
            distances = Distances(self._perfect_structure,
                                  site.frac_coords,
                                  self.dist_tol)
            result.append(distances.atom_idx_at_center(specie=site.specie))
        return result

    @property
    def removed_indices(self):
        result = []
        for p, d in enumerate(self.p_to_d):
            try:
                if self.d_to_p[d] != p:
                    result.append(p)
            except (IndexError, TypeError):
                result.append(p)
        return sorted(result)

    @property
    def inserted_indices(self):
        result = []
        for d, p in enumerate(self.d_to_p):
            try:
                if self.p_to_d[p] != d:
                    result.append(d)
            except (IndexError, TypeError):
                result.append(d)
        return sorted(result)

    @property
    def defect_center_coord(self):
        coords = []
        for v in self.removed_indices:
            coords.append(self._perfect_structure[v].frac_coords)
        for i in self.inserted_indices:
            coords.append(self._defect_structure[i].frac_coords)

        lattice = self._perfect_structure.lattice
        repr_coords = coords[0]
        translated_coords = [list(repr_coords)]
        for c in coords[1:]:
            _, trans = lattice.get_distance_and_image(repr_coords, c)
            translated_coords.append([c[i] + trans[i] for i in range(3)])
        return np.average(translated_coords, axis=0) % 1

    def neighboring_atom_indices(self, cutoff_factor=None):
        distances = []
        for v in self.removed_indices:
            distances.append(Distances(self._defect_structure,
                                       self._perfect_structure[v].frac_coords,
                                       self.dist_tol))
        for i in self.inserted_indices:
            distances.append(Distances(self._defect_structure,
                                       self._defect_structure[i].frac_coords,
                                       self.dist_tol))
        result = set()
        for d in distances:
            result.update(d.coordination(cutoff_factor=cutoff_factor)
                          .neighboring_atom_indices)
        return sorted(list(result))

    def make_site_diff(self):
        removed = {}
        for idx in self.removed_indices:
            site = self._perfect_structure[idx]
            removed[idx] = site.species_string, tuple(site.frac_coords)

        inserted = {}
        for idx in self.inserted_indices:
            site = self._defect_structure[idx]
            inserted[idx] = site.species_string, tuple(site.frac_coords)

        return SiteDiff(
            removed=removed, inserted=inserted, mapping=self.atom_mapping)


@dataclass
class SiteDiff(MSONable):
    removed: Dict[int, Tuple[str, Tuple[float, float, float]]]
    inserted: Dict[int, Tuple[str, Tuple[float, float, float]]]
    mapping: Dict[int, int]


#
# def make_defect_type(vacancies, interstitials, lattice, same_dist_criterion):
#     if len(vacancies) == 1 and len(interstitials) == 0:
#         return DefectType.vacancy
#     if len(vacancies) == 0 and len(interstitials) == 1:
#         return DefectType.interstitial
#
#
# class DefectType(MSONable, ExtendedEnum):
#     vacancy = "vacancy"
#     interstitial = "interstitial"
#     substituted = "substituted"
#     complex = "complex"

