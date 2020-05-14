# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass
from pymatgen import Structure, Element
from monty.json import MSONable


@dataclass
class Site(MSONable):
    wyckoff_letter: str
    site_symmetry: str
    equivalent_atoms: List[int]


@dataclass
class SupercellInfo(MSONable):
    structure: Structure
    transform_matrix: List[List[int]]
    sites: Dict[str, Site]



