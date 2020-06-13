# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import List, Tuple

import numpy as np
from pydefect.defaults import defaults
from pydefect.input_maker.defect import SimpleDefect
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.input_maker.defect_set import DefectSet
from pydefect.input_maker.supercell_info import SupercellInfo
from pymatgen import Structure, IStructure


class DefectEntriesMaker:
    def __init__(self, supercell_info: SupercellInfo, defect_set: DefectSet):
        self.supercell_info = supercell_info
        self.defect_entries = set()

        for defect in defect_set:
            structure, perturbed_structure, coords, site_symmetry = \
                self._create_defect_structures(defect)

            for charge in defect.charges:
                self.defect_entries.add(
                    DefectEntry(name=defect.name,
                                charge=charge,
                                structure=structure,
                                perturbed_structure=perturbed_structure,
                                site_symmetry=site_symmetry,
                                defect_center=coords))

    def _create_defect_structures(self,
                                  defect: SimpleDefect,
                                  ) -> Tuple[IStructure, IStructure,
                                             Tuple[float, float, float], str]:

        structure = copy_to_structure(self.supercell_info.structure)

        if defect.out_atom[0] is "i":
            index = int(defect.out_atom[1:]) - 1
            site = self.supercell_info.interstitials[index]
            cutoff = self.supercell_info.interstitial_coords(index).cutoff
            coords = self.supercell_info.interstitials[index].frac_coords
        else:
            site = self.supercell_info.sites[defect.out_atom]
            cutoff = self.supercell_info.coords(defect.out_atom).cutoff
            removed_site_index = site.equivalent_atoms[0]
            coords = structure.pop(removed_site_index).frac_coords

        perturbed_structure = perturb_structure(structure, coords, cutoff)

        if defect.in_atom:
            add_atom_to_structure(structure, defect.in_atom, coords)
            add_atom_to_structure(perturbed_structure, defect.in_atom, coords)

        return to_istructure(structure), to_istructure(perturbed_structure), \
               tuple(coords), site.site_symmetry


def copy_to_structure(structure: IStructure) -> Structure:
    return Structure.from_dict(structure.as_dict())


def to_istructure(structure: Structure) -> IStructure:
    return IStructure.from_dict(structure.as_dict())


def add_atom_to_structure(structure: Structure, elem: str, coords: List[float]):
    """In-place atom insertion to structure
    https://stackoverflow.com/questions/2361426/get-the-first-item-from-an-iterable-that-matches-a-condition
    """
    try:
        idx = next(i for i, s in enumerate(structure) if str(s.specie) == elem)
    except StopIteration:
        idx = len(structure)

    structure.insert(idx, elem, coords)


def perturb_structure(structure: Structure, center: List[float], cutoff: float
                      ) -> Structure:
    """ structure perturbation

    Args:
        structure: pmg Structure class object
        center: Fractional coordinates of a central position.
        cutoff: Radius of a sphere in which atoms are perturbed.
    """
    result = structure.copy()
    cartesian_coords = structure.lattice.get_cartesian_coords(center)
    neighboring_atoms = structure.get_sites_in_sphere(
        pt=cartesian_coords, r=cutoff, include_index=True)

#    assert cutoff < min(structure.lattice.lengths) / 2.0

    # neighboring_atom is composed of (PeriodicSite, distance, index)
    for _, _, site_index in neighboring_atoms:
        vector = random_3d_vector(defaults.displace_distance)
        result.translate_sites(site_index, vector, frac_coords=False)

    return result


def random_3d_vector(distance):
    """Random 3d vector with uniform spherical distribution with 0 <= norm <= 1.
    stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    """
    phi = np.random.uniform(0, np.pi * 2)
    cos_theta = np.random.uniform(-1, 1)
    theta = np.arccos(cos_theta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    return np.array([x, y, z]) * np.random.uniform(high=distance)
