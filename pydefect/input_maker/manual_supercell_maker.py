# -*- coding: utf-8 -*-
#  Copyright (c) 2023. Distributed under the terms of the MIT License.

from typing import List, Dict

import yaml
from numpy.linalg import det
from pydefect.defaults import defaults
from pydefect.input_maker.supercell import Supercell, Supercells
from pydefect.input_maker.supercell_info import SupercellInfo, SimpleSite
from pymatgen.core import IStructure
from vise.util.logger import get_logger
from vise.util.structure_symmetrizer import StructureSymmetrizer

logger = get_logger(__name__)


class ManualSupercellMaker:

    def __init__(self,
                 input_structure: IStructure,
                 sites: Dict[str, SimpleSite],
                 matrix: List[List[int]] = None,
                 symprec: float = defaults.symmetry_length_tolerance,
                 angle_tolerance: float = defaults.symmetry_angle_tolerance,
                 **supercell_kwargs):

        symmetrizer = StructureSymmetrizer(input_structure,
                                           symprec=symprec,
                                           angle_tolerance=angle_tolerance)
        if input_structure != symmetrizer.primitive:
            logger.warning(
                "The input structure is different from the primitive.")

        if matrix:
            self.supercell = Supercell(input_structure, matrix)
        else:
            supercells = Supercells(input_structure, **supercell_kwargs)
            self.supercell = supercells.most_isotropic_supercell

        for site in sites.values():
            multiplicity = int(round(det(self.supercell.matrix)))
            site.site_index *= multiplicity

        self.supercell_info = SupercellInfo(
            structure=self.supercell.structure,
            space_group=symmetrizer.spglib_sym_data["international"],
            transformation_matrix=self.supercell.matrix,
            sites=sites,
            unitcell_structure=input_structure)


def make_sites_from_yaml_file(filename: str):
    try:
        with open(filename) as file:
            d = yaml.safe_load(file)
    except TypeError:
        logger.warning("Yaml file format for creating site info is incorrect.")
        raise TypeError
    return {k: SimpleSite(**v) for k, v in d.items()}
