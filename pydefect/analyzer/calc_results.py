# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List

import numpy as np
from pymatgen import IStructure


@dataclass
class CalcResults:
    structure: IStructure
    site_symmetry: str
    energy: float
    magnetization: float
    eigenvalues: np.array
    kpoint_coords: List[List[float]]
    kpoint_weights: List[float]
    electrostatic_potential: List[float]
    vbm: float
    cbm: float
    volume: float
    fermi_level: float


