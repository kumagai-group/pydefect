# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import Union, List

import numpy as np
from monty.json import MSONable
from pymatgen import IStructure, Structure
from vise.util.structure_symmetrizer import StructureSymmetrizer

from pydefect.defaults import defaults


@dataclass
class DefectStructure(MSONable):
    name: str
    charge: int
    initial_structure: IStructure
    final_structure: IStructure
    initial_symmetry: str
    final_symmetry: str
    defect_center: Union[int, List[float]]

    @property
    def anchor_atom_index(self):
        """ Returns an index of atom that is the farthest from the defect.

        This atom is assumed not to displace in the defective supercell, and
        so used for analyzing local structure around a defect.
        Note that only the first occurrence is returned when using argmax.
        docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.argmax.html
        """
        lattice = self.initial_structure.lattice
        all_coords = self.initial_structure.frac_coords
        dist_set = lattice.get_all_distances(self.defect_center, all_coords)[0]
        return int(np.argmax(dist_set))




def symmetrize_defect_structure(structure: IStructure,
                                anchor_atom_idx: int,
                                anchor_atom_coord) -> Structure:
    ss = StructureSymmetrizer(structure, defaults.symmetry_length_tolerance)
    result = ss.primitive.copy()
    offset = result[anchor_atom_idx].frac_coords - anchor_atom_coord
    result.translate_sites(list(range(len(structure))), offset)
    return result

