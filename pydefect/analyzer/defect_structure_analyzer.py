# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
from pymatgen import IStructure, Structure
from vise.util.structure_symmetrizer import StructureSymmetrizer

from pydefect.defaults import defaults
from pydefect.util.structure_tools import Distances


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

        d_to_p = []
        for site in defective_structure:
            distances = Distances(perfect_structure, site.frac_coords)
            d_to_p.append(distances.mapped_atom_idx(specie=site.specie))

        p_to_d = []
        for site in perfect_structure:
            distances = Distances(defective_structure, site.frac_coords)
            p_to_d.append(distances.mapped_atom_idx(specie=site.specie))

        inserted_indices = []
        for d, p in enumerate(d_to_p):
            try:
                if p_to_d[p] != d:
                    inserted_indices.append(d)
            except (IndexError, TypeError):
                inserted_indices.append(d)
        self.inserted_indices = sorted(inserted_indices)

        vacancy_indices = []
        for p, d in enumerate(p_to_d):
            try:
                if d_to_p[d] != p:
                    vacancy_indices.append(p)
            except (IndexError, TypeError):
                vacancy_indices.append(p)
        self.vacancy_indices = sorted(vacancy_indices)

        self.atom_mapping = \
            {d: p for d, p in enumerate(d_to_p) if d not in inserted_indices}

    @property
    def defect_center_coord(self):
        coords = []
        for v in self.vacancy_indices:
            coords.append(self._perfect_structure[v].frac_coords)
        for i in self.inserted_indices:
            coords.append(self._defective_structure[i].frac_coords)
        return np.average(coords, axis=0)


def symmetrize_defect_structure(structure: IStructure,
                                anchor_atom_idx: int,
                                anchor_atom_coord) -> Structure:
    ss = StructureSymmetrizer(structure, defaults.symmetry_length_tolerance)
    result = ss.primitive.copy()
    offset = result[anchor_atom_idx].frac_coords - anchor_atom_coord
    result.translate_sites(list(range(len(structure))), offset)
    return result

