# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass

import numpy as np
from monty.json import MSONable
from pymatgen import Lattice


@dataclass
class Ewald(MSONable):
    """Container class for anisotropic Ewald parameter and related properties

    lattice (Lattice):
        The given lattice.
    dielectric_tensor (3x3 np.array):
        Static dielectric tensor where the directions are compatible
        with the lattice.
    ewald_param (float):
        A parameter used for evaluating Ewald sum.
    accuracy (float):
        product of cutoff radius of G-vector and gaussian FWHM. Increasing this
        value, calculation will be more accurate, but slower.
    """
    lattice: Lattice
    dielectric_tensor: np.ndarray
    ewald_param: float
    accuracy: float

    @property
    def volume(self):
        return self.lattice.volume

    @property
    def rec_lattice_matrix(self):
        return self.lattice.reciprocal_lattice.matrix


