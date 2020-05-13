# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from typing import List
import numpy as np
from numpy.linalg import det

from pymatgen import Structure


class Supercell:
    def __init__(self, input_structure: Structure, matrix: List[List[int]]):
        self.matrix = matrix
        self.structure = input_structure * matrix

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
                 input_structure: Structure,
                 min_num_atoms: int,
                 max_num_atoms: int):
        self._input_structure = input_structure
        min_det = min_num_atoms / len(input_structure) - 1e-5
        max_det = max_num_atoms / len(input_structure) + 1e-5
        matrix = np.eye(3, dtype=int)

        self.supercells = []

        while True:
            if det(matrix) >= min_det:
                self.supercells.append(Supercell(input_structure, matrix))
            if det(matrix) > max_det:
                break
            idx = self.min_length_index(matrix)
            matrix[idx][idx] += 1

    def min_length_index(self, matrix):
        abc = (self._input_structure * matrix).lattice.lengths
        return abc.index(min(abc))

    @property
    def isotropic_supercell(self):
        return min(self.supercells, key=lambda s: s.isotropy)



