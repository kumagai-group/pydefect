# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from typing import Optional, List

import numpy as np
from numpy.linalg import det
from pydefect.defaults import defaults
from pydefect.input_maker.supercell import Supercell, TetragonalSupercells, \
    Supercells
from pydefect.input_maker.supercell_info import SupercellInfo
from pydefect.util.error_classes import NotPrimitiveError, SupercellError
from pymatgen.core import IStructure
from vise.util.centering import Centering
from vise.util.logger import get_logger
from vise.util.structure_symmetrizer import StructureSymmetrizer, Site

logger = get_logger(__name__)


class SupercellMaker:
    def __init__(self,
                 primitive_structure: IStructure,
                 matrix_to_conv_cell: Optional[List[List[int]]] = None,
                 symprec: float = defaults.symmetry_length_tolerance,
                 angle_tolerance: float = defaults.symmetry_angle_tolerance,
                 raise_error: bool = True,
                 **supercell_kwargs):

        self.primitive_structure = primitive_structure
        self.symmetrizer = StructureSymmetrizer(primitive_structure,
                                           symprec=symprec,
                                           angle_tolerance=angle_tolerance)
        if primitive_structure != self.symmetrizer.primitive:
            logger.warning(
                "The input structure is different from the primitive one,"
                "which might be due to the difference of symprec used in"
                "the pydefect and unitcell conversion.")
            logger.warning("\n".join([
                "Input lattice:",
                f"{primitive_structure.lattice}", "",
                "Primitive structure lattice:",
                f"{self.symmetrizer.primitive.lattice}", "",
                "Input structure:",
                f"{primitive_structure}", "",
                "Primitive structure:",
                f"{self.symmetrizer.primitive}"]))
            if raise_error:
                raise NotPrimitiveError

        self.sg_symbol = self.symmetrizer.spglib_sym_data["international"]
        self.conv_structure = self.symmetrizer.conventional
        crystal_system, center = str(self.symmetrizer.bravais)

        centering = Centering(center)
        self.conv_multiplicity = centering.conv_multiplicity
        self.conv_trans_mat = centering.primitive_to_conv

        self._matrix = matrix_to_conv_cell
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
        multiplicity = int(round(det(self.transformation_matrix)))

        sites = {}
        for site_name, site in self.symmetrizer.sites.items():
            atoms = []
            for i in site.equivalent_atoms:
                atoms.extend(list(range(i*multiplicity, (i+1)*multiplicity)))
            sites[site_name] = Site(element=site.element,
                                    wyckoff_letter=site.wyckoff_letter,
                                    site_symmetry=site.site_symmetry,
                                    equivalent_atoms=atoms)

        self.supercell_info = \
            SupercellInfo(structure=self.supercell.structure,
                          space_group=self.sg_symbol,
                          transformation_matrix=self.transformation_matrix,
                          sites=sites,
                          unitcell_structure=self.primitive_structure)

    @property
    def transformation_matrix(self) -> List[List[int]]:
        """
        Definition of transformation_matrix is different from that in spglib.
        https://spglib.github.io/spglib/definition.html#transformation-matrix-boldsymbol-p-and-origin-shift-boldsymbol-p
        This is because we need to keep the consistency with the pymatgen.
        See also the unittest.
        """
        matrix = np.dot(self.supercell.matrix, self.conv_trans_mat).astype(int)
        return matrix.tolist()

