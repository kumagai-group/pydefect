# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np
from typing import List, Tuple
from pymatgen import Structure, IStructure

from pydefect.input_maker.defect import SimpleDefect
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.input_maker.defect_set import DefectSet
from pydefect.input_maker.supercell_info import SupercellInfo, Site

from pydefect.defaults import defaults


class DefectEntriesMaker:
    def __init__(self, supercell_info: SupercellInfo, defect_set: DefectSet):
        self.supercell_info = supercell_info
        self.defect_entries = set()

        for defect in defect_set:
            site = supercell_info.sites[defect.out_atom]
            _, cutoff = supercell_info.coords(defect.out_atom)
            structure, perturbed_structure = \
                self._create_defect_structure(defect, site, cutoff)

            for charge in defect.charges:
                self.defect_entries.add(DefectEntry(
                    name=defect.name,
                    charge=charge,
                    initial_structure=structure,
                    perturbed_initial_structure=perturbed_structure,
                    initial_site_symmetry=site.site_symmetry))

    def _create_defect_structure(self,
                                 defect: SimpleDefect,
                                 site: Site,
                                 cutoff: float
                                 ) -> Tuple[IStructure, IStructure]:
        removed_site_index = site.equivalent_atoms[0]
        structure = copy_to_structure(self.supercell_info.structure)
        coords = structure.pop(removed_site_index).frac_coords
        perturbed_structure = perturb_structure(structure, coords, cutoff)
        if defect.in_atom:
            add_atom_to_structure(structure, defect.in_atom, coords)
            add_atom_to_structure(perturbed_structure, defect.in_atom, coords)

        return to_istructure(structure), to_istructure(perturbed_structure)


def copy_to_structure(structure: IStructure) -> Structure:
    return Structure.from_dict(structure.as_dict())


def to_istructure(structure: Structure) -> IStructure:
    return IStructure.from_dict(structure.as_dict())


def add_atom_to_structure(structure: Structure,
                          elem: str,
                          coords: List[float]) -> None:
    """In-place atom insertion to structure"""
    index = len(structure)
    for i, site in enumerate(structure):
        if str(site.specie) == elem:
            index = i
            break

    structure.insert(index, elem, coords)


def perturb_structure(structure: Structure,
                      center: List[float],
                      cutoff: float) -> Structure:
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

    assert cutoff < min(structure.lattice.lengths) / 2.0

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
