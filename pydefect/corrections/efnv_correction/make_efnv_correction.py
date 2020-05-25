# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
from numpy import dot, cross
from numpy.linalg import norm

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_structure_analyzer import DefectStructureAnalyzer
from pydefect.corrections.efnv_correction.efnv_correction import \
    ExtendedFnvCorrection, DefectSite
from pydefect.corrections.efnv_correction.ewald import Ewald
from pydefect.input_maker.defect_entry import DefectEntry


def make_efnv_correction(defect_entry: DefectEntry,
                         calc_results: CalcResults,
                         perfect_calc_results: CalcResults,
                         dielectric_tensor: np.array):
    """
    Notes:
    (1) The formula written in YK2014 are divided by 4pi to keep the SI unit.
    (2) When assuming an element charge locate at the defect_coords and
        angstrom for length, relative dielectric tensor, Multiply
        elementary_charge * 1e10 / epsilon_0 = 180.95128169876497
        to make potential in V.
    """

    structure_analyzer = DefectStructureAnalyzer(
        calc_results.structure, perfect_calc_results.structure)
    defect_coord = structure_analyzer.defect_center_coord
    lattice = calc_results.structure.lattice
    ewald = Ewald(lattice, dielectric_tensor)
    charge = defect_entry.charge
    point_charge_correction = - ewald.lattice_energy * charge ** 2

    sites = []
    for d, p in enumerate(structure_analyzer.atom_mapping):
        specie = calc_results.structure[d].specie
        distance = structure_analyzer.distance_from_center(d)
        potential = (calc_results.electrostatic_potential[d] -
                     perfect_calc_results.electrostatic_potential[p])

        coord = calc_results.structure[d].frac_coords
        rel_coord = [x - y for x, y in zip(coord, defect_coord)]
        pc_potential = ewald.atomic_site_potential(rel_coord) * charge

        sites.append(DefectSite(specie, distance, potential, pc_potential))

    return ExtendedFnvCorrection(
        charge=charge,
        point_charge_correction=point_charge_correction,
        defect_region_radius=calc_max_sphere_radius(lattice),
        sites=sites)


def calc_max_sphere_radius(lattice) -> float:
    """Calculate Maximum radius of a sphere fitting inside the unit cell.

    Calculate three distances between two parallel planes using
    (a_i x a_j) . a_k / |a_i . a_j| """
    distances = np.zeros(3, dtype=float)
    for i in range(3):
        a_i_a_j = cross(lattice[i - 2], lattice[i - 1])
        a_k = lattice[i]
        distances[i] = abs(dot(a_i_a_j, a_k)) / norm(a_i_a_j)
    return max(distances) / 2.0
