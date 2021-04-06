# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Optional

import numpy as np
import yaml
from monty.json import MSONable
from pydefect.analyzer.defect_structure_comparator import \
    DefectStructureComparator
from pymatgen import Site
from pymatgen.core import IStructure
from vise.util.mix_in import ToJsonFileMixIn
from vise.util.structure_symmetrizer import StructureSymmetrizer
from vise.util.typing import Coords


@dataclass(frozen=True)
class PerturbedSite(MSONable):
    element: str
    distance: float
    initial_coords: Coords
    perturbed_coords: Coords
    displacement: float

    def __str__(self):
        return


@dataclass(frozen=True)
class DefectEntry(MSONable, ToJsonFileMixIn):
    name: str
    charge: int
    structure: IStructure
    site_symmetry: str
    defect_center: Coords
    perturbed_structure: Optional[IStructure] = None
    perturbed_sites: Optional[Tuple[PerturbedSite, ...]] = None
    perturbed_site_symmetry: Optional[str] = None

    @property
    def anchor_atom_index(self) -> int:
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

    def __str__(self):
        center = ", ".join([f"{c:6.3f}" for c in self.defect_center])
        perturbed = []
        for i in self.perturbed_site_indices:
            i_site: Site = self.structure[i]
            p_site = self.perturbed_structure[i]
            perturbed.append([i_site.specie, ])

        return f""" -- defect entry info
name: {self.full_name}
site symmetry: {self.site_symmetry}
defect center: ({center})    
perturbed sites:
elem dist   initial_coords             perturbed_coords
  He 5.00 ( 0.000,  0.000,  0.500) -> ( 0.000,  0.000,  0.501) 
  He 5.00 ( 0.000,  0.500,  0.000) -> ( 0.000,  0.501,  0.000) 
  He 5.00 ( 0.500,  0.000,  0.000) -> ( 0.501,  0.000,  0.000)"""


def make_defect_entry(name: str,
                      charge: int,
                      perfect_structure: IStructure,
                      defect_structure: IStructure):

    analyzer = DefectStructureComparator(perfect_structure, defect_structure)

    species = []
    frac_coords = []
    for d, p in enumerate(analyzer.p_to_d):
        if p is None:
            site = defect_structure[d]
        else:
            site = perfect_structure[p]
        species.append(site.specie)
        frac_coords.append(site.frac_coords)

    initial_structure = IStructure(perfect_structure.lattice,
                                   species, frac_coords)
    symmetrizer = StructureSymmetrizer(initial_structure)

    return DefectEntry(name=name,
                       charge=charge,
                       structure=initial_structure,
                       site_symmetry=symmetrizer.point_group,
                       defect_center=tuple(analyzer.defect_center_coord))
