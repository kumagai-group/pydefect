# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass

from pymatgen import IStructure, Structure
from vise.util.structure_symmetrizer import StructureSymmetrizer

from pydefect.defaults import defaults


@dataclass
class DefectStructure:
    initial_structure: IStructure
    final_structure: IStructure
    initial_symmetry: str
    final_symmetry: str
    defect_center: int


def symmetrize_defect_structure(structure: IStructure,
                                anchor_atom_idx: int,
                                anchor_atom_coord) -> Structure:
    ss = StructureSymmetrizer(structure, defaults.symmetry_length_tolerance)
    result = ss.primitive.copy()
    offset = result[anchor_atom_idx].frac_coords - anchor_atom_coord
    result.translate_sites(list(range(len(structure))), offset)
    return result

