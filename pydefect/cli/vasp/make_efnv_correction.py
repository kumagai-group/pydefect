# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import Optional, Tuple

import numpy as np

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_structure_comparator import \
    DefectStructureComparator
from pydefect.corrections.defect_region import HalfMinFaceDistanceDefectRegion, DefectRegion
from pydefect.corrections.efnv_correction import \
    ExtendedFnvCorrection, PotentialSite
from pydefect.corrections.ewald import Ewald
from pydefect.defaults import defaults
from pydefect.util.error_classes import SupercellError, \
    NoCalculatedPotentialSiteError


Coords = Tuple[float, float, float]


def make_efnv_correction(charge: float,
                         calc_results: CalcResults,
                         perfect_calc_results: CalcResults,
                         dielectric_tensor: np.array,
                         defect_coords: Optional[Coords] = None,
                         accuracy: float = defaults.ewald_accuracy,
                         defect_region: DefectRegion = None,
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
    if defect_region is None:
        defect_region = HalfMinFaceDistanceDefectRegion(sample_radius_ratio=1.0)

    sites, rel_coords, defect_coords = \
        make_sites(calc_results, perfect_calc_results, defect_coords)

    lattice = calc_results.structure.lattice

    radius = defect_region.defect_region_radius(lattice.matrix)

    ewald = Ewald(lattice.matrix, dielectric_tensor, accuracy=accuracy)
    point_charge_correction = \
        0.0 if not charge else - ewald.lattice_energy * charge ** 2

    has_calculated_sites = False
    for site, rel_coord in zip(sites, rel_coords):
        if calc_all_sites is True or site.distance > radius:
            has_calculated_sites = True
            if charge == 0:
                site.pc_potential = 0
            else:
                site.pc_potential = (ewald.atomic_site_potential(rel_coord)
                                     * charge * unit_conversion)

    if has_calculated_sites is False:
        raise NoCalculatedPotentialSiteError(
            "Change the spherical radius of defect region. "
            f"Now {radius:4.2f}Å is set.")

    return ExtendedFnvCorrection(
        charge=charge,
        point_charge_correction=point_charge_correction * unit_conversion,
        defect_region_radius=radius,
        defect_region=defect_region,
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



