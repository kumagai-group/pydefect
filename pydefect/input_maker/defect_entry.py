# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Optional

import numpy as np
import yaml
from monty.json import MSONable, MontyDecoder
from pydefect.analyzer.defect_structure_comparator import \
    DefectStructureComparator
from pydefect.util.coords import pretty_coords
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

    @classmethod
    def from_dict(cls, d):
        return cls(element=d["element"],
                   distance=d["distance"],
                   initial_coords=tuple(d["initial_coords"]),
                   perturbed_coords=tuple(d["perturbed_coords"]),
                   displacement=d["displacement"])

    def __str__(self):
        elem = f"{self.element:>4}"
        dist = f"{self.distance:4.2f}"
        ini_coords = pretty_coords(self.initial_coords)
        p_coords = pretty_coords(self.perturbed_coords)
        disp = f"{self.displacement:7.2f}"
        return f"{elem} {dist} {ini_coords} -> {p_coords} {disp}"


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

    @classmethod
    def from_dict(cls, d):
        d = {k: MontyDecoder().process_decoded(v) for k, v in d.items() if not k.startswith("@")}
        if d["perturbed_sites"]:
            d["perturbed_sites"] = tuple([MontyDecoder().process_decoded(d) for d in d["perturbed_sites"]])
        d["defect_center"] = tuple(d["defect_center"])
        return cls(**d)

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
        result = f""" -- defect entry info
name: {self.full_name}
site symmetry: {self.site_symmetry}
defect center: ({", ".join([f"{c:6.3f}" for c in self.defect_center])})
perturbed sites:
elem dist   initial_coords             perturbed_coords         displacement
"""
        if self.perturbed_sites:
            result += "\n".join([i.__str__() for i in self.perturbed_sites])

        return result


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
