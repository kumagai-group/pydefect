# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from dataclasses import dataclass
from monty.json import MSONable
from pymatgen import IStructure


@dataclass(frozen=True)
class DefectEntry(MSONable):
    name: str
    charge: int
    structure: IStructure
    perturbed_structure: IStructure
    site_symmetry: str

    @property
    def perturbed_site_indices(self):
        result = []
        for i, (s1, s2) in enumerate(zip(self.structure,
                                         self.perturbed_structure)):
            if (s1.frac_coords != s2.frac_coords).any():
                result.append(i)
        return result
