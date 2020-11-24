# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import Optional

import numpy as np
from numpy.linalg import inv
from pydefect.util.structure_tools import Distances
from pymatgen import IStructure, Structure
from vise.util.logger import get_logger
from vise.util.structure_symmetrizer import StructureSymmetrizer

logger = get_logger(__name__)


class DefectStructureAnalyzer:
    def __init__(self,
                 defective_structure: IStructure,
                 perfect_structure: IStructure):
        """
        Atoms in the final structure are shifted such that the farthest atom
        from the defect is placed at the same place with that in the perfect
        supercell.
        """
        self._defective_structure = defective_structure
        self._perfect_structure = perfect_structure
        self.p_to_d = self.make_p_to_d()
        self.d_to_p = self.make_d_to_p()

    @property
    def atom_mapping(self):
        return {d: p for d, p in enumerate(self.d_to_p)
                if d not in self.inserted_indices}

    def make_p_to_d(self):
        result = []
        for site in self._perfect_structure:
            distances = Distances(self._defective_structure, site.frac_coords)
            result.append(distances.atom_idx_at_center(specie=site.specie))
        return result

    def make_d_to_p(self):
        result = []
        for site in self._defective_structure:
            distances = Distances(self._perfect_structure, site.frac_coords)
            result.append(distances.atom_idx_at_center(specie=site.specie))
        return result

    @property
    def vacancy_indices(self):
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
        for v in self.vacancy_indices:
            coords.append(self._perfect_structure[v].frac_coords)
        for i in self.inserted_indices:
            coords.append(self._defective_structure[i].frac_coords)

        lattice = self._perfect_structure.lattice
        repr_coords = coords[0]
        translated_coords = [list(repr_coords)]
        for c in coords[1:]:
            _, trans = lattice.get_distance_and_image(repr_coords, c)
            translated_coords.append([c[i] + trans[i] for i in range(3)])
        return np.average(translated_coords, axis=0) % 1

    def distance_from_center(self, index):
        frac_coords = self._defective_structure[index].frac_coords
        result, _ = self._defective_structure.lattice.get_distance_and_image(
            self.defect_center_coord, frac_coords)
        return result

    @property
    def neighboring_atom_indices(self):
        result = set()
        for v in self.vacancy_indices:
            distances = Distances(self._defective_structure,
                                  self._perfect_structure[v].frac_coords)
            result.update(distances.coordination().neighboring_atom_indices)
        for i in self.inserted_indices:
            distances = Distances(self._defective_structure,
                                  self._defective_structure[i].frac_coords)
            result.update(distances.coordination().neighboring_atom_indices)

        return sorted(list(result))

    def displacements(self, anchor_atom_idx):
        anchor_fcoords_after = self._defective_structure[anchor_atom_idx]
        anchor_fcoords_before = self._perfect_structure[self.d_to_p[anchor_atom_idx]]
        offset = (anchor_fcoords_before.frac_coords
                  - anchor_fcoords_after.frac_coords)
        result = []
        lattice = self._defective_structure.lattice
        for d, p in enumerate(self.d_to_p):
            if p is None:
                result.append(None)
            else:
                displacement = (self._perfect_structure[p].frac_coords
                                - self._defective_structure[d].frac_coords
                                - offset)
                displacement = fold_frac_coords(displacement)
                result.append(lattice.get_cartesian_coords(displacement))
        return result


def fold_frac_coords(frac_coords):
    return [(frac_coord + 0.5) % 1 - 0.5 for frac_coord in frac_coords]


def symmetrize_defect_structure(structure_symmetrizer: StructureSymmetrizer,
                                anchor_atom_idx: Optional[int] = None,
                                anchor_atom_coord: Optional[np.ndarray] = None
                                ) -> Structure:
    result = structure_symmetrizer.structure.copy()

    logger.info(f"The symmetry is {structure_symmetrizer.point_group}")

    origin_shift = structure_symmetrizer.spglib_sym_data["origin_shift"]
    inv_trans_mat = inv(structure_symmetrizer.spglib_sym_data["transformation_matrix"])
    coords = structure_symmetrizer.spglib_sym_data["std_positions"]
    new_coords = []
    for i in range(len(result)):
        new_coords.append(np.dot(inv_trans_mat, (coords[i] - origin_shift)))

    if anchor_atom_idx:
        offset = new_coords[anchor_atom_idx] - anchor_atom_coord
    else:
        offset = np.array([0.0]*3)
    new_coords = np.array(new_coords) - offset
    for i, coords in zip(result, new_coords):
        i.frac_coords = coords % 1
    return result

