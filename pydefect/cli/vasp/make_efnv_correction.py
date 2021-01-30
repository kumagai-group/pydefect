# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

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


def make_efnv_correction(charge: int,
                         calc_results: CalcResults,
                         perfect_calc_results: CalcResults,
                         dielectric_tensor: np.array,
                         accuracy: float = defaults.ewald_accuracy,
                         unit_conversion: float = 180.95128169876497):
    """
    Notes:
    (1) The formula written in YK2014 need to be divided by 4pi in the SI unit.
    (2) When assuming an element charge locate at the defect_coords and
        angstrom for length, relative dielectric tensor, Multiply
        elementary_charge * 1e10 / epsilon_0 = 180.95128169876497
        to make potential in V.
    """
    assert calc_results.structure.lattice == perfect_calc_results.structure.lattice

    structure_analyzer = DefectStructureComparator(
        calc_results.structure, perfect_calc_results.structure)
    defect_coords = structure_analyzer.defect_center_coord
    lattice = calc_results.structure.lattice
    ewald = Ewald(lattice.matrix, dielectric_tensor, accuracy=accuracy)
    point_charge_correction = \
        0.0 if not charge else - ewald.lattice_energy * charge ** 2

    defect_region_radius = calc_max_sphere_radius(lattice.matrix)

    sites = []
    for d, p in structure_analyzer.atom_mapping.items():
        specie = str(calc_results.structure[d].specie)
        frac_coords = calc_results.structure[d].frac_coords
        distance, _ = lattice.get_distance_and_image(defect_coords, frac_coords)
        pot = calc_results.potentials[d] - perfect_calc_results.potentials[p]

        coord = calc_results.structure[d].frac_coords
        rel_coord = [x - y for x, y in zip(coord, defect_coords)]
        if distance <= defect_region_radius:
            pc_potential = None
        else:
            if charge == 0:
                pc_potential = 0
            else:
                pc_potential = ewald.atomic_site_potential(rel_coord) * charge
                pc_potential *= unit_conversion

        sites.append(PotentialSite(specie, distance, pot, pc_potential))

    return ExtendedFnvCorrection(
        charge=charge,
        point_charge_correction=point_charge_correction * unit_conversion,
        defect_region_radius=defect_region_radius,
        sites=sites,
        defect_coords=tuple(defect_coords))


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
