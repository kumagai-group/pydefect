# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np

from pydefect.corrections.efnv_correction.ewald import Ewald


class AnisotropicPotential:
    """
    Notes:
    (1) The formula written in YK2014 are divided by 4pi to keep the SI unit.

    (2) When assuming an element charge locate at the defect_coords and
        angstrom for length, relative dielectric tensor, Multiply
        elementary_charge * 1e10 / epsilon_0 = 180.95128169876497
        to make potential in V.
    """

    def __init__(self,
                 defect_coord: np.ndarray,
                 target_coords: np.ndarray,
                 ewald: Ewald):
        model_pot = []
        for coord in target_coords:
            shift = [x - y for x, y in zip(coord, defect_coord)]
            real_part = ewald.ewald_real(include_self=True, shift=shift)
            rec_part = ewald.ewald_rec(shift)
            model_pot.append(real_part + rec_part + ewald.diff_pot)

        self.model_pot = model_pot

        real_part = ewald.ewald_real(include_self=False, shift=[0, 0, 0])
        rec_part = ewald.ewald_rec([0, 0, 0])
        self.lattice_energy = (real_part + rec_part + ewald.diff_pot + ewald.self_pot) / 2
