# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass
from pymatgen import Structure, Element
from monty.json import MSONable
from numpy.linalg import det
from pydefect.util.structure_tools import Distances


@dataclass
class Site(MSONable):
    wyckoff_letter: str
    site_symmetry: str
    equivalent_atoms: List[int]


@dataclass
class SupercellInfo(MSONable):
    structure: Structure
    space_group: str
    transform_matrix: List[List[int]]
    sites: Dict[str, Site]

    def coords(self, name):
        site = self.sites[name]
        coord = list(self.structure[site.equivalent_atoms[0]].frac_coords)
        distances = Distances(self.structure, coord)
        return distances.coord_distances_and_cutoff()

