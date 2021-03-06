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