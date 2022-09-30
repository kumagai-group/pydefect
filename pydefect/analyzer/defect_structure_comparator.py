# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
from monty.json import MSONable
from pydefect.defaults import defaults
from pydefect.util.structure_tools import Distances
from pymatgen.core import IStructure, Structure
from vise.util.typing import Coords


class DefectStructureComparator:
    def __init__(self,
                 defect_structure: IStructure,
                 perfect_structure: IStructure,
                 dist_tol: float = defaults.dist_tol):
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

    def _atom_projection(self, structure_from, structure_to, specie=True):
        result = []
        for site in structure_from:
            distances = Distances(structure_to,
                                  site.frac_coords,
                                  self.dist_tol)
            specie = site.specie if specie else None
            result.append(distances.atom_idx_at_center(specie=specie))
        return result

    def make_p_to_d(self):
        return self._atom_projection(
            self._perfect_structure, self._defect_structure)

    def make_d_to_p(self):
        return self._atom_projection(
            self._defect_structure, self._perfect_structure)

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
        removed_sites = [self._perfect_structure[x]
                         for x in self.removed_indices]
        inserted_sites = [self._defect_structure[x]
                          for x in self.inserted_indices]

        try:
            removed_str = Structure.from_sites(removed_sites)
        except ValueError:
            removed_str = None

        try:
            inserted_str = Structure.from_sites(inserted_sites)
        except ValueError:
            inserted_str = None

        if inserted_str and removed_str:
            r_to_i = self._atom_projection(removed_str, inserted_str, specie=False)
            i_to_r = self._atom_projection(inserted_str, removed_str, specie=False)
        else:
            r_to_i = [None]*len(removed_sites)
            i_to_r = [None]*len(inserted_sites)

        mapping, removed_mapping, inserted_mapping = [], [], []
        for x, y in enumerate(r_to_i):
            if y and x == i_to_r[y]:
                mapping.append((self.removed_indices[x], self.inserted_indices[y]))
                removed_mapping.append(self.removed_indices[x])
                inserted_mapping.append(self.inserted_indices[y])

        removed, removed_by_sub = [], []
        for idx in self.removed_indices:
            site = self._perfect_structure[idx]
            frac_coords = tuple([float(fc) for fc in site.frac_coords])
            val = idx, site.species_string, frac_coords
            if idx in removed_mapping:
                removed_by_sub.append(val)
            else:
                removed.append(val)

        inserted, inserted_by_sub = [], []
        for idx in self.inserted_indices:
            site = self._defect_structure[idx]
            frac_coords = tuple([float(fc) for fc in site.frac_coords])
            val = idx, site.species_string, frac_coords
            if idx in inserted_mapping:
                inserted_by_sub.append(val)
            else:
                inserted.append(val)

        return SiteDiff(removed=removed,
                        inserted=inserted,
                        removed_by_sub=removed_by_sub,
                        inserted_by_sub=inserted_by_sub)


SiteInfo = Tuple[int, str, Coords]


@dataclass
class SiteDiff(MSONable):
    removed: List[SiteInfo]
    inserted: List[SiteInfo]
    removed_by_sub: List[SiteInfo]
    inserted_by_sub: List[SiteInfo]

    @classmethod
    def from_dict(cls, d):
        result = super().from_dict(d)
        removed, inserted, removed_by_sub, inserted_by_sub = [], [], [], []

        for i in result.removed:
            removed.append((i[0], i[1], tuple(i[2])))
        for i in result.inserted:
            inserted.append((i[0], i[1], tuple(i[2])))
        for i in result.removed_by_sub:
            removed_by_sub.append((i[0], i[1], tuple(i[2])))
        for i in result.inserted_by_sub:
            inserted_by_sub.append((i[0], i[1], tuple(i[2])))

        return cls(removed, inserted, removed_by_sub, inserted_by_sub)

    @property
    def is_complex_defect(self):
        return (len(self.removed) + len(self.inserted)
                + len(self.removed_by_sub)) != 1

    @property
    def is_vacancy(self):
        return (len(self.removed) == 1 and len(self.inserted) == 0 and
                len(self.removed_by_sub) == 0 and len(self.inserted_by_sub) == 0)

    @property
    def is_interstitial(self):
        return (len(self.removed) == 0 and len(self.inserted) == 1 and
                len(self.removed_by_sub) == 0 and len(self.inserted_by_sub) == 0)

    @property
    def is_substituted(self):
        return (len(self.removed) == 0 and len(self.inserted) == 0 and
                len(self.removed_by_sub) == 1 and len(self.inserted_by_sub) == 1)

    @property
    def is_no_diff(self):
        return not (self.removed or self.inserted
                    or self.removed_by_sub or self.inserted_by_sub)

