# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from dataclasses import dataclass
from monty.json import MSONable
from pymatgen import IStructure


@dataclass(frozen=True)
class DefectEntry(MSONable):
    name: str
    charge: int
    initial_structure: IStructure
    perturbed_initial_structure: IStructure
    initial_site_symmetry: str

    @property
    def perturbed_site_indices(self):
        result = []
        for i, (s1, s2) in enumerate(zip(self.initial_structure,
                                         self.perturbed_initial_structure)):
            if (s1.frac_coords != s2.frac_coords).any():
                result.append(i)
        return result
