# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import List

import numpy as np


def sanitize_matrix(matrix: List[int]) -> List[List[int]]:

    if len(matrix) == 9:
        return [matrix[:3], matrix[3:6], matrix[6:]]
    elif len(matrix) == 3:
        result = np.eye(3, dtype=int)
        for i in range(3):
            result[i, i] = matrix[i]
        return result.tolist()
    elif len(matrix) == 1:
        result = np.eye(3, dtype=int)
        for i in range(3):
            result[i, i] = matrix[0]
        return result.tolist()
    else:
        raise ValueError(f"Matrix element length {len(matrix)} is improper.")


def str_int_to_int(x):
    try:
        return int(x)
    except ValueError:
        return x
