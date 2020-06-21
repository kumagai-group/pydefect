# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np
import yaml
from monty.json import MSONable
from pydefect.util.mix_in import ToJsonFileMixIn
from pymatgen import IStructure


@dataclass(frozen=True)
class DefectEntry(MSONable, ToJsonFileMixIn):
    name: str
    charge: int
    structure: IStructure
    perturbed_structure: IStructure
    site_symmetry: str
    defect_center: Tuple[float, float, float]

    @property
    def anchor_atom_index(self):
        """ Returns an index of atom that is the farthest from the defect.

         Only the first occurrence is returned when using argmax.
         docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.argmax.html
         """
        lattice = self.structure.lattice
        all_coords = self.structure.frac_coords
        dist_set = lattice.get_all_distances(self.defect_center, all_coords)[0]
        return int(np.argmax(dist_set))

    @property
    def anchor_atom_coords(self) -> np.ndarray:
        return self.structure[self.anchor_atom_index].frac_coords

    @property
    def perturbed_site_indices(self):
        result = []
        for i, (s1, s2) in enumerate(zip(self.structure,
                                         self.perturbed_structure)):
            if (s1.frac_coords != s2.frac_coords).any():
                result.append(i)
        return result

    @property
    def full_name(self):
        return "_".join([self.name, str(self.charge)])

    def to_prior_info(self, filename):
        d = {"charge": self.charge}
        Path(filename).write_text(yaml.dump(d))

