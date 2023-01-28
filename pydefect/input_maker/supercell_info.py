# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Union

from monty.json import MSONable
from numpy.linalg import det
from pydefect.database.database import electronegativity, oxidation_state
from pydefect.util.structure_tools import Distances
from pymatgen.core import IStructure
from vise.util.logger import get_logger
from vise.util.mix_in import ToJsonFileMixIn
from vise.util.structure_symmetrizer import Site

logger = get_logger(__name__)


@dataclass(frozen=True)
class Interstitial(MSONable):
    frac_coords: List[float]
    site_symmetry: str = None
    info: str = None


@dataclass
class SimpleSite(MSONable):
    element: str
    site_index: int
    site_symmetry: str = None
    wyckoff_letter: str = None

    @property
    def pprint_equiv_atoms(self):
        return None

    @property
    def equivalent_atoms(self):
        return [self.site_index]


@dataclass(frozen=True)
class SupercellInfo(MSONable, ToJsonFileMixIn):
    structure: IStructure
    space_group: str
    transformation_matrix: List[List[int]]
    sites: Dict[str, Union[Site, SimpleSite]]
    interstitials: List[Interstitial] = field(default_factory=list)
    unitcell_structure: IStructure = None

    @classmethod
    def from_dict(cls, d):
        if "transform_matrix" in d:
            d["transformation_matrix"] = d.pop("transform_matrix")
        return super().from_dict(d)

    def coords(self, name):
        site = self.sites[name]
        site_idx = site.equivalent_atoms[0]
        coord = list(self.structure[site_idx].frac_coords)
        distances = Distances(self.structure, coord)
        return distances.coordination()

    def interstitial_coords(self, idx: int):
        interstitial = self.interstitials[idx]
        distances = Distances(self.structure, interstitial.frac_coords)
        return distances.coordination(include_on_site=True)

    @property
    def multiplicity(self):
        return int(round(det(self.transformation_matrix)))

    @property
    def _transform_matrix_str(self):
        return '  '.join(str(x) for x in self.transformation_matrix)

    def _frac_coords(self, site):
        repr_atom_idx = site.equivalent_atoms[0]
        frac_coords = self.structure[repr_atom_idx].frac_coords
        return "{0[0]:9.7f}  {0[1]:9.7f}  {0[2]:9.7f}".format(frac_coords)

    def __str__(self):
        lines = [f"Space group: {self.space_group}",
                 f"Transformation matrix: {self._transform_matrix_str}",
                 f"Cell multiplicity: {self.multiplicity}", ""]

        for name, site in self.sites.items():
            elem = stripe_numbers(name)
            coordination = self.coords(name).distance_dict
            lines.append(f"   Irreducible element: {name}")
            lines.append(f"        Wyckoff letter: {site.wyckoff_letter}")
            lines.append(f"         Site symmetry: {site.site_symmetry}")
            lines.append(f"         Cutoff radius: {self.coords(name).cutoff}")
            lines.append(f"          Coordination: {coordination}")
            lines.append(f"      Equivalent atoms: {site.pprint_equiv_atoms}")
            lines.append(f"Fractional coordinates: {self._frac_coords(site)}")
            lines.append(f"     Electronegativity: {electronegativity(elem)}")
            lines.append(f"       Oxidation state: {oxidation_state(elem)}")
            lines.append("")
        lines.append("-- interstitials")

        for idx, site in enumerate(self.interstitials):
            coordination = self.interstitial_coords(idx).distance_dict

            lines.append(f"#{idx + 1}")
            if site.info:
                lines.append(f"                  Info: {site.info}")

            frac = \
                "{0[0]:9.7f}  {0[1]:9.7f}  {0[2]:9.7f}".format(site.frac_coords)
            lines.append(f"Fractional coordinates: {frac}")
            lines.append(f"         Site symmetry: {site.site_symmetry}")
            lines.append(f"          Coordination: {coordination}")
            lines.append("")

        return "\n".join(lines)


def stripe_numbers(name):
    return ''.join([i for i in name if i.isalpha()])
