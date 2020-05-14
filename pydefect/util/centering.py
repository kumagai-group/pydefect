# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np

from vise.util.enum import ExtendedEnum


class Centering(ExtendedEnum):
    P = "P"
    A = "A"
    C = "C"
    R = "R"
    I = "I"
    F = "F"

    @property
    def conv_to_primitive(self):
        """ Note that the matrices are transposed from those written in spglib
        web sites due to the difference of matrix convention.

        See unittests for sure.

        https://spglib.github.io/spglib/definition.html
        """
        if self is self.P:
            return np.eye(3)
        elif self is self.A:
            return np.array([[1,    0,    0],
                             [0,  1/2,  1/2],
                             [0, -1/2,  1/2]])
        elif self is self.C:
            return np.array([[1/2, -1/2, 0],
                             [1/2,  1/2, 0],
                             [  0,    0, 1]])
        elif self is self.R:
            return np.array([[ 2/3,  1/3,  1/3],
                             [-1/3,  1/3,  1/3],
                             [-1/3, -2/3,  1/3]])
        elif self is self.I:
            return np.array([[-1/2,  1/2,  1/2],
                             [ 1/2, -1/2,  1/2],
                             [ 1/2,  1/2, -1/2]])
        elif self is self.F:
            return np.array([[  0,  1/2, 1/2],
                             [1/2,    0, 1/2],
                             [1/2,  1/2,   0]])
        else:
            raise NotImplementedError(f"Define matrix for {self}")

    @property
    def primitive_to_conv(self):
        return np.linalg.inv(self.conv_to_primitive)

    @property
    def conv_multiplicity(self):
        return int(round(1 / np.linalg.det(self.conv_to_primitive)))

    # return matrix

