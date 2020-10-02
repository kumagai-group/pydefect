# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from typing import Optional, List

import numpy as np
from pydefect.input_maker.supercell import Supercell, TetragonalSupercells, \
    Supercells
from pydefect.input_maker.supercell_info import SupercellInfo
from pydefect.util.error_classes import NotPrimitiveError, SupercellError
from pymatgen import IStructure
from vise.util.centering import Centering
from vise.util.logger import get_logger
from vise.util.structure_symmetrizer import StructureSymmetrizer

logger = get_logger(__name__)


class SupercellMaker:
    def __init__(self,
                 primitive_structure: IStructure,
                 matrix: Optional[List[List[int]]] = None,
                 **supercell_kwargs):

        self.primitive_structure = primitive_structure
        symmetrizer = StructureSymmetrizer(primitive_structure)
        if primitive_structure != symmetrizer.primitive:
            logger.warning(f"Input structure: {primitive_structure}")
            logger.warning(f"Primitive structure:{symmetrizer.primitive}")
            raise NotPrimitiveError

        self.sg = symmetrizer.sg_number
        self.sg_symbol = symmetrizer.spglib_sym_data["international"]
        self.conv_structure = symmetrizer.conventional
        crystal_system, center = str(symmetrizer.bravais)

        centering = Centering(center)
        self.conv_multiplicity = centering.conv_multiplicity
        self.conv_trans_mat = centering.primitive_to_conv

        self._matrix = matrix
        self._supercell_kwargs = supercell_kwargs

        self._generate_supercell(crystal_system)
        self._generate_supercell_info()

    def _generate_supercell(self, crystal_system):
        if self._matrix:
            self.supercell = Supercell(self.conv_structure, self._matrix)
        else:
            if crystal_system == "t":
                supercells = TetragonalSupercells(self.conv_structure,
                                                  **self._supercell_kwargs)
            else:
                supercells = Supercells(self.conv_structure,
                                        **self._supercell_kwargs)

            self.supercell = supercells.most_isotropic_supercell

    def _generate_supercell_info(self):
        symmetrizer = StructureSymmetrizer(self.supercell.structure)
        if symmetrizer.sg_number != self.sg:
            raise SupercellError

        self.supercell_info = SupercellInfo(self.supercell.structure,
                                            self.sg_symbol,
                                            self.transformation_matrix,
                                            symmetrizer.sites,
                                            unitcell_structure=self.primitive_structure)

    @property
    def transformation_matrix(self):
        matrix = np.dot(self.conv_trans_mat, self.supercell.matrix).astype(int)
        return matrix.tolist()

