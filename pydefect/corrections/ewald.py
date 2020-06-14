# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from functools import reduce
from math import sqrt, pow, ceil
from typing import List, Optional

import numpy as np
from monty.json import MSONable
from numpy import sqrt, dot, pi, exp, cos
from numpy.linalg import norm, inv, det
from pydefect.defaults import defaults
from scipy.special import erfc
from scipy.stats import mstats


def grid_number(lattice_vectors: np.ndarray, max_length: float):
    a = [ceil(max_length / norm(lattice_vectors[i])) for i in range(3)]
    return a[0] * a[1] * a[2]


class Ewald(MSONable):
    """Container class for anisotropic Ewald parameter and related properties

    ewald_param:
        A parameter used for evaluating Ewald sum.
    accuracy:
        product of cutoff radius of G-vector and gaussian FWHM. Increasing this
        value, calculation will be more accurate, but slower.
    """

    def __init__(self,
                 lattice: np.ndarray,
                 dielectric_tensor: np.ndarray,
                 accuracy: float = defaults.ewald_accuracy,
                 ewald_param: Optional[float] = None):
        self.lattice = lattice
        self.rec_lattice = inv(self.lattice).T * 2 * pi
        self.volume = det(self.lattice)
        self.cube_root_vol = pow(self.volume, 1 / 3)

        self.dielectric_tensor = dielectric_tensor
        self.det_epsilon = det(dielectric_tensor)
        self.root_epsilon = sqrt(self.det_epsilon)
        self.epsilon_inv = inv(dielectric_tensor)

        self.accuracy = accuracy
        if ewald_param:
            self.ewald_param = ewald_param
        else:
            l_r = mstats.gmean([norm(v) for v in self.lattice])
            l_g = mstats.gmean([norm(v) for v in self.rec_lattice])
            self.ewald_param = sqrt(l_g / l_r / 2) * self.cube_root_vol / self.root_epsilon

        # Modified Ewald parameter which is the gamma in the Eqs in YK2014.
        self.mod_ewald_param = self.ewald_param / self.cube_root_vol * self.root_epsilon
        # 2nd term in Eq.(14) in YK2014, caused by finite gaussian charge.
        self.diff_pot = -0.25 / self.volume / self.mod_ewald_param ** 2

    def atomic_site_potential(self, rel_coord):
        real_part = self.ewald_real(include_self=True, shift=rel_coord)
        rec_part = self.ewald_rec(rel_coord)
        return real_part + rec_part + self.diff_pot

    @property
    def lattice_energy(self):
        real_part = self.ewald_real(include_self=False, shift=[0, 0, 0])
        rec_part = self.ewald_rec([0, 0, 0])
        return (real_part + rec_part + self.diff_pot + self.self_pot) / 2

    @property
    def self_pot(self):
        return - self.mod_ewald_param / (2.0 * pi * sqrt(pi * self.det_epsilon))

    def ewald_real(self, include_self: bool, shift: List[float]) -> float:
        summed = 0
        for r in self.r_lattice_set(include_self, shift):
            root_r_inv_epsilon_r = sqrt(reduce(dot, [r.T, self.epsilon_inv, r]))
            summed += erfc(self.mod_ewald_param * root_r_inv_epsilon_r) / root_r_inv_epsilon_r
        return summed / (4 * pi * self.root_epsilon)

    def ewald_rec(self, coord) -> float:
        cart_coord = dot(coord, self.lattice)
        summed = 0
        for g in self.g_lattice_set():
            g_epsilon_g = reduce(dot, [g.T, self.dielectric_tensor, g])
            summed += (exp(- g_epsilon_g / 4 / self.mod_ewald_param ** 2)
                       / g_epsilon_g * cos(dot(g, cart_coord)))
        return summed / self.volume

    def r_lattice_set(self,
                      include_self: bool = True,
                      shift: List[float] = None) -> np.ndarray:
        xyz = self.xyz(self.r_vector_nums, shift)
        if not include_self:
            xyz = np.delete(xyz, int((len(xyz) - 1) / 2), 0)
        return dot(xyz, self.lattice)

    def g_lattice_set(self) -> np.ndarray:
        xyz = self.xyz(self.g_vector_nums)
        xyz = np.delete(xyz, int((len(xyz) - 1) / 2), 0)
        return dot(xyz, self.rec_lattice)

    @staticmethod
    def xyz(nums: List[int], shift: List[float] = None) -> np.ndarray:
        shift = shift or [0, 0, 0]
        x = np.arange(-nums[0], nums[0] + 1, 1) - shift[0]
        y = np.arange(-nums[1], nums[1] + 1, 1) - shift[1]
        z = np.arange(-nums[2], nums[2] + 1, 1) - shift[2]
        return np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)

    @property
    def r_vector_nums(self) -> List[int]:
        max_length = self.accuracy / self.mod_ewald_param
        return [ceil(max_length / norm(self.lattice[i])) for i in range(3)]

    @property
    def g_vector_nums(self) -> List[int]:
        max_length = 2 * self.mod_ewald_param * self.accuracy
        return [ceil(max_length / norm(self.rec_lattice[i])) for i in range(3)]

