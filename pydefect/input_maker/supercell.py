# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from copy import deepcopy
from functools import reduce
from itertools import product
from typing import List, Tuple

import numpy as np
from numpy.linalg import det
from pydefect.util.error_classes import SupercellError
from pymatgen.core import IStructure


class Supercell:
    def __init__(self, input_structure: IStructure, matrix: List[List[int]]):
        self.matrix = matrix
        self.structure = input_structure * matrix
        self.lattice = self.structure.lattice

    @property
    def isotropy(self):
        lengths = self.structure.lattice.lengths
        average = np.average(lengths)
        return sum([abs(length - average) for length in lengths]) / 3 / average

    @property
    def average_angle(self):
        return sum(self.structure.lattice.angles) / 3


class Supercells:
    def __init__(self,
                 input_structure: IStructure,
                 min_num_atoms: int = 50,
                 max_num_atoms: int = 250):
        self.input_structure = input_structure
        min_det = min_num_atoms / len(input_structure) - 1e-5
        max_det = max_num_atoms / len(input_structure) + 1e-5

        self.supercells = []
        matrix = np.eye(3, dtype=int)
        for i in range(50):
            if det(matrix) > max_det:
                break
            if det(matrix) > min_det:
                self.supercells.append(Supercell(input_structure,
                                                 matrix.tolist()))
            matrix = self.incremented_matrix(matrix)

    def incremented_matrix(self, matrix: np.ndarray):
        abc = (self.input_structure * matrix).lattice.lengths
        min_indices = [i for i, val in enumerate(abc) if val == min(abc)]
        new_matrix = deepcopy(matrix)
        for index in min_indices:
            new_matrix[index][index] += 1
        return new_matrix

    @property
    def most_isotropic_supercell(self):
        try:
            return min(self.supercells, key=lambda s: s.isotropy)
        except ValueError:
            raise SupercellError


class RhombohedralSupercells(Supercells):
    def incremented_matrix(self, matrix: np.ndarray):
        alpha = (self.input_structure * matrix).lattice.angles[0]

        fan_out = np.array([[1,  1, -1], [-1,  1,  1], [1, -1,  1]])
        close_in = np.array([[1, 0, 1], [1, 1, 0], [0, 1, 1]])

        if alpha < 90:
            return np.dot(matrix, fan_out)
        else:
            return np.dot(matrix, close_in)


class TetragonalSupercells(Supercells):
    def incremented_matrix(self, matrix: np.ndarray):
        a, _, c = (self.input_structure * matrix).lattice.lengths
        new_matrix = deepcopy(matrix)
        if a > c:
            new_matrix[2][2] += 1
            return new_matrix
        else:
            ab_submatrix_det = round(det(matrix[0:2, 0:2]))
            x, y = self.next_x_y_combination(ab_submatrix_det)
            new_matrix[0:2, 0:2] = self.matrix_from_x_y(x, y)
            return new_matrix

    @staticmethod
    def matrix_from_x_y(x: int, y: int) -> np.ndarray:
        """ """
        assert x >= 1
        assert y >= 0

        expand_matrix = np.array([[x, 0], [0, x]])
        if y == 0:
            rot_matrix = np.eye(2)
        elif y == 1:
            rot_matrix = np.array([[1, 1], [-1, 1]])
        else:
            matrix_list = [np.array([[1, 1], [-1, 1]])] * y
            rot_matrix = reduce(np.dot, matrix_list)

        return np.dot(expand_matrix, rot_matrix)

    @staticmethod
    def next_x_y_combination(num: int) -> Tuple[int, int]:
        """ Find x, y with minimum of x ** 2 * 2^y larger than num."""
        pairs = {(x, y): x * x * 2 ** y
                 for x, y in product(range(1, 10), range(8))}
        filtered = {k: v for k, v in pairs.items() if v > num}
        return next(iter(sorted(filtered, key=lambda u: pairs[u])))


