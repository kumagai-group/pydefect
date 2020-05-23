# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from math import sqrt, pow, ceil

import numpy as np
from monty.json import MSONable
from numpy import dot, cross, pi
from numpy.linalg import norm, inv, det
from pymatgen import IStructure
from scipy.stats import mstats


def calc_max_sphere_radius(lattice_vectors: np.ndarray) -> float:
    """ Calculate three distances between two parallel planes using
        (a_i x a_j) . a_k / |a_i . a_j|
    """
    distances = np.zeros(3, dtype=float)
    for i in range(3):
        a_i_a_j = cross(lattice_vectors[i - 2], lattice_vectors[i - 1])
        a_k = lattice_vectors[i]
        distances[i] = abs(dot(a_i_a_j, a_k)) / norm(a_i_a_j)
    return max(distances) / 2.0


def grid_number(lattice_vectors: np.ndarray, max_length: float):
    a = [ceil(max_length / norm(lattice_vectors[i])) for i in range(3)]
    return a[0] * a[1] * a[2]


@dataclass
class Ewald(MSONable):
    """Container class for anisotropic Ewald parameter and related properties

    ewald_param:
        A parameter used for evaluating Ewald sum.
    accuracy:
        product of cutoff radius of G-vector and gaussian FWHM. Increasing this
        value, calculation will be more accurate, but slower.
    """
    lattice: np.ndarray
    dielectric_tensor: np.ndarray
    ewald_param: float
    accuracy: float

    @property
    def reciprocal_lattice(self):
        return inv(self.lattice).T * 2 * pi

    @property
    def volume(self):
        return det(self.lattice)


def make_ewald(structure: IStructure,
               dielectric_tensor: np.ndarray,
               convergence: float = 1.05,
               accuracy: float = 25.0):
    """Get optimized ewald parameter.

    determine initial ewald parameter to satisfy following:
    max_int(Real) = max_int(Reciprocal)
    in neighbor_lattices function.

    Left term:
    max_int(Real) = 2 * x * Y  / l_r, where x, Y, and l_r are ewald,
    prod_cutoff_fwhm, and axis length of real lattice, respectively.

    Right term:
    max_int(reciprocal) = Y  / (x * l_g)
    where l_g is axis length of reciprocal lattice, respectively.
    Then, x = sqrt(l_g / l_r / 2)
    gmean : geometric mean,  (a1 * a2 * a3)^(1/3)

    initial_ewald_param (float):
        Initial guess of parameter.
    convergence (float):
        If 1/convergence < n_(real)/n_(reciprocal) < convergence,
        where n_(real) and n_(reciprocal) is number of real lattices
        and reciprocal lattices, finishes optimization and
        returns ewald_param.
    accuracy (float):
        product of cutoff radius of G-vector and gaussian FWHM.
        Increasing this value, calculation will be more accurate, but slower.
    """
    root_det_dielectric = np.sqrt(np.linalg.det(dielectric_tensor))
    lattice = structure.lattice.matrix
    reciprocal_lattice = structure.lattice.reciprocal_lattice.matrix
    cube_root_vol = pow(structure.lattice.volume, 1 / 3)

    l_r = mstats.gmean([norm(v) for v in lattice])
    l_g = mstats.gmean([norm(v) for v in reciprocal_lattice])
    ewald_param = sqrt(l_g / l_r / 2) * cube_root_vol / root_det_dielectric

    for i in range(10):
        ewald = ewald_param / cube_root_vol * root_det_dielectric

        max_r_vector_norm = accuracy / ewald
        real_num = grid_number(lattice, max_r_vector_norm)

        max_g_vector_norm = 2 * ewald * accuracy
        rec_num = grid_number(reciprocal_lattice, max_g_vector_norm)

        if 1 / convergence < real_num / rec_num < convergence:
            return Ewald(lattice, dielectric_tensor, ewald_param, accuracy)
        else:
            ewald_param *= (real_num / rec_num) ** (1 / 6)
    else:
        raise ValueError("The initial ewald param may not be adequate.")



