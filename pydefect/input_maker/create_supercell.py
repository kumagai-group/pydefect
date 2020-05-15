# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from typing import Optional, List
import numpy as np
from pymatgen import Structure
from vise.util.structure_symmetrizer import StructureSymmetrizer
from itertools import groupby
from collections import defaultdict

from pydefect.input_maker.supercell import Supercell, TetragonalSupercells, \
    Supercells
from pydefect.util.error_classes import NotPrimitiveError, SupercellError
from pydefect.util.centering import Centering
from pydefect.input_maker.supercell_info import SupercellInfo, Site


class CreateSupercell:
    def __init__(self,
                 input_structure: Structure,
                 matrix: Optional[List[List[int]]] = None,
                 **supercell_kwargs):
        symmetrizer = StructureSymmetrizer(input_structure)
        self.sg = symmetrizer.sg_number
        self.sg_symbol = symmetrizer.spglib_sym_data["international"]
        if input_structure.lattice != symmetrizer.primitive.lattice:
            raise NotPrimitiveError

        self.conv_structure = symmetrizer.conventional
        crystal_system, center = str(symmetrizer.bravais)
        centering = Centering(center)
        self.conv_multiplicity = centering.conv_multiplicity
        self.conv_trans_mat = centering.primitive_to_conv

        self._generate_supercell(crystal_system, matrix, supercell_kwargs)
        self._generate_supercell_info()

    def _generate_supercell(self, crystal_system, matrix, supercell_kwargs):
        if matrix:
            self.supercell = Supercell(self.conv_structure, matrix)
        else:
            if crystal_system == "t":
                supercells = TetragonalSupercells(self.conv_structure,
                                                  **supercell_kwargs)
            else:
                supercells = Supercells(self.conv_structure, **supercell_kwargs)

            self.supercell = supercells.most_isotropic_supercell

    def _generate_supercell_info(self):
        symmetrizer = StructureSymmetrizer(self.supercell.structure)
        if symmetrizer.sg_number != self.sg:
            raise SupercellError

        wyckoffs = symmetrizer.spglib_sym_data["wyckoffs"]
        equivalent_atoms = symmetrizer.spglib_sym_data["equivalent_atoms"]
        site_symmetries = symmetrizer.spglib_sym_data["site_symmetry_symbols"]

        idx_equiv = sorted(enumerate(equivalent_atoms), key=lambda x: x[1])
        sites = {}
        element_idx_dict = defaultdict(int)

        for _, equiv_sites in groupby(idx_equiv, lambda x: x[1]):
            equiv_site_list = list(equiv_sites)
            repr_idx = equiv_site_list[0][0]
            element = self.supercell.structure[repr_idx].specie.name
            element_idx_dict[element] += 1
            index = str(element_idx_dict[str(element)])
            name = element + index
            sites[name] = Site(element=element,
                               wyckoff_letter=wyckoffs[repr_idx],
                               site_symmetry=site_symmetries[repr_idx],
                               equivalent_atoms=[s[0] for s in equiv_site_list])

        self.supercell_info = SupercellInfo(self.supercell.structure,
                                            self.sg_symbol,
                                            self.transform_matrix,
                                            sites)

    @property
    def transform_matrix(self):
        return np.dot(self.conv_trans_mat, self.supercell.matrix).tolist()