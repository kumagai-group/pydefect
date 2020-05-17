# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from dataclasses import dataclass

from monty.json import MSONable
from pymatgen import IStructure


@dataclass
class DefectEntry(MSONable):
    name: str
    charge: int
    initial_structure: IStructure
    perturbed_initial_structure: IStructure
    initial_site_symmetry: str



