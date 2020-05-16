# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from dataclasses import dataclass

from monty.json import MSONable
from pymatgen import Structure


@dataclass
class DefectEntry(MSONable):
    name: str
    charge: int
    initial_structure: Structure
    perturbed_initial_structure: Structure
    initial_site_symmetry: str



