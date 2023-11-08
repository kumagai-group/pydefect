# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import Optional

import numpy as np
from numpy import dot, cross
from numpy.linalg import norm

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_structure_comparator import \
    DefectStructureComparator
from pydefect.corrections.efnv_correction import \
    ExtendedFnvCorrection, PotentialSite
from pydefect.corrections.ewald import Ewald
from pydefect.defaults import defaults
from pydefect.util.error_classes import SupercellError, \
    NoCalculatedPotentialSiteError
from vise.util.logger import get_logger
from vise.util.typing import Coords

logger = get_logger(__name__)


def make_efnv_correction(charge: float,
                         calc_results: CalcResults,
                         perfect_calc_results: CalcResults,
                         dielectric_tensor: np.array,
                         defect_coords: Optional[Coords] = None,
                         accuracy: float = defaults.ewald_accuracy,
                         defect_region_radius: float = None,
                         calc_all_sites: bool = False,
                         unit_conversion: float = 180.95128169876497):
    """
    Notes:
    (1) The formula written in YK2014 need to be divided by 4pi in the SI unit.
    (2) When assuming an element charge locate at the defect_coords and
        angstrom for length, relative dielectric tensor, Multiply
        elementary_charge * 1e10 / epsilon_0 = 180.95128169876497
        to make potential in V.
    """
    sites, rel_coords, defect_coords = \
        make_sites(calc_results, perfect_calc_results, defect_coords)

    lattice = calc_results.structure.lattice
    ewald = Ewald(lattice.matrix, dielectric_tensor, accuracy=accuracy)
    point_charge_correction = \
        0.0 if not charge else - ewald.lattice_energy * charge ** 2
    if defect_region_radius is None:
        defect_region_radius = calc_max_sphere_radius(lattice.matrix)

    has_calculated_sites = False
    for site, rel_coord in zip(sites, rel_coords):
        if calc_all_sites is True or site.distance > defect_region_radius:
            has_calculated_sites = True
            if charge == 0:
                site.pc_potential = 0
            else:
                site.pc_potential = (ewald.atomic_site_potential(rel_coord)
                                     * charge * unit_conversion)

    if has_calculated_sites is False:
        raise NoCalculatedPotentialSiteError(
            "Change the spherical radius of defect region manually. "
            f"Now {defect_region_radius:4.2f}Å is set.")

    return ExtendedFnvCorrection(
        charge=charge,
        point_charge_correction=point_charge_correction * unit_conversion,
        defect_region_radius=defect_region_radius,
        sites=sites,
        defect_coords=tuple(defect_coords))


def make_sites(calc_results, perfect_calc_results, defect_coords):
    if calc_results.structure.lattice != perfect_calc_results.structure.lattice:
        raise SupercellError("The lattice constants for defect and perfect "
                             "models are different")
    structure_analyzer = DefectStructureComparator(
        calc_results.structure, perfect_calc_results.structure)
    if defect_coords is None:
        defect_coords = structure_analyzer.defect_center_coord
    lattice = calc_results.structure.lattice
    sites, rel_coords = [], []

    for d, p in structure_analyzer.atom_mapping.items():
        specie = str(calc_results.structure[d].specie)
        frac_coords = calc_results.structure[d].frac_coords
        distance, _ = lattice.get_distance_and_image(defect_coords, frac_coords)
        pot = calc_results.potentials[d] - perfect_calc_results.potentials[p]
        sites.append(PotentialSite(specie, distance, pot, None))
        coord = calc_results.structure[d].frac_coords
        rel_coords.append([x - y for x, y in zip(coord, defect_coords)])

    return sites, rel_coords, defect_coords


def calc_max_sphere_radius(lattice_matrix) -> float:
    """Calculate Maximum radius of a sphere fitting inside the unit cell.

    Calculate three distances between two parallel planes using
    (a_i x a_j) . a_k / |a_i . a_j| """
    distances = np.zeros(3, dtype=float)
    for i in range(3):
        a_i_a_j = cross(lattice_matrix[i - 2], lattice_matrix[i - 1])
        a_k = lattice_matrix[i]
        distances[i] = abs(dot(a_i_a_j, a_k)) / norm(a_i_a_j)
    return max(distances) / 2.0
