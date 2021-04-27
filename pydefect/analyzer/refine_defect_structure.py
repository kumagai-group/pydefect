# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import numpy as np
from numpy.linalg import inv
from pydefect.analyzer.defect_structure_info import logger
from pydefect.defaults import defaults
from pymatgen.core import Structure
from vise.util.structure_symmetrizer import StructureSymmetrizer


def refine_defect_structure(structure: Structure,
                            anchor_atom_index: int = None,
                            anchor_atom_coords: np.ndarray = None):
    symmetrizer = StructureSymmetrizer(structure,
                                       defaults.symmetry_length_tolerance,
                                       defaults.symmetry_angle_tolerance)
    result = structure.copy()
    spglib_data = symmetrizer.spglib_sym_data

    origin_shift = spglib_data["origin_shift"]
    inv_trans_mat = inv(spglib_data["transformation_matrix"])
    coords = spglib_data["std_positions"]
    # This transformation is key for refinement.
    new_coords = np.array([np.dot(inv_trans_mat, (coords[i] - origin_shift))
                           for i in range(len(result))])
    if anchor_atom_index:
        offset = new_coords[anchor_atom_index] - anchor_atom_coords
        new_coords = new_coords - offset
    for i, coords in zip(result, new_coords):
        i.frac_coords = coords % 1

    if result != structure:
        logger.info(f"Structure is refined to the one with point group "
                    f"{symmetrizer.point_group}.")
        return result
    else:
        logger.info(f"Refined structure is same as before, so do nothing.")
