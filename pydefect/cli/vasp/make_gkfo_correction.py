# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np

from pydefect.analyzer.calc_results import CalcResults
from pydefect.corrections.efnv_correction import \
    ExtendedFnvCorrection, PotentialSite
from pydefect.corrections.ewald import Ewald
from pydefect.corrections.gkfo_correction import GkfoCorrection
from pydefect.defaults import defaults


def make_gkfo_correction(efnv_correction: ExtendedFnvCorrection,
                         additional_charge: int,
                         final_calc_results: CalcResults,
                         initial_calc_results: CalcResults,
                         diele_tensor: np.array,
                         ion_clamped_diele_tensor: np.array,
                         accuracy: float = defaults.ewald_accuracy,
                         unit_conversion: float = 180.95128169876497
                         ) -> GkfoCorrection:
    assert final_calc_results.structure == initial_calc_results.structure
    assert abs(additional_charge) == 1

    defect_coords = efnv_correction.defect_coords
    lattice = initial_calc_results.structure.lattice
    ewald_ele = Ewald(lattice.matrix, ion_clamped_diele_tensor,
                      accuracy=accuracy)
    defect_region_radius = efnv_correction.defect_region_radius

    pc_2nd_term = - ewald_ele.lattice_energy

    gkfo_sites = []
    for index, site in enumerate(initial_calc_results.structure):
        specie = site.specie
        dist, _ = lattice.get_distance_and_image(site.frac_coords, defect_coords)
        pot = (final_calc_results.potentials[index]
               - initial_calc_results.potentials[index])

        rel_coord = [x - y for x, y in zip(site.frac_coords, defect_coords)]
        if dist <= defect_region_radius:
            pc_potential = None
        else:
            pc_potential = (ewald_ele.atomic_site_potential(rel_coord)
                            * additional_charge * unit_conversion)

        gkfo_sites.append(PotentialSite(specie, dist, pot, pc_potential))

    return GkfoCorrection(
        init_efnv_correction=efnv_correction,
        additional_charge=additional_charge,
        pc_2nd_term=pc_2nd_term * unit_conversion,
        gkfo_sites=gkfo_sites,
        ave_dielectric_tensor=np.trace(diele_tensor) / 3,
        ave_electronic_dielectric_tensor=np.trace(ion_clamped_diele_tensor) / 3)



