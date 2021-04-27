# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from typing import Optional

import numpy as np
from numpy.linalg import inv
from pydefect.analyzer.defect_structure_info import logger
from pymatgen.core import Structure
from vise.util.structure_symmetrizer import StructureSymmetrizer


def symmetrize_defect_structure(structure_symmetrizer: StructureSymmetrizer,
                                anchor_atom_idx: Optional[int] = None,
                                anchor_atom_coord: Optional[np.ndarray] = None
                                ) -> Structure:
    result = structure_symmetrizer.structure.copy()
    spglib_data = structure_symmetrizer.spglib_sym_data

    origin_shift = spglib_data["origin_shift"]
    inv_trans_mat = inv(spglib_data["transformation_matrix"])
    coords = spglib_data["std_positions"]

    new_coords = np.array([np.dot(inv_trans_mat, (coords[i] - origin_shift))
                           for i in range(len(result))])
    if anchor_atom_idx:
        offset = new_coords[anchor_atom_idx] - anchor_atom_coord
        new_coords = new_coords - offset

    for i, coords in zip(result, new_coords):
        i.frac_coords = coords % 1
    return result
