# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Optional

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
    vbm_info: float
    cbm_info: float
    fermi_level: float
    electronic_conv: Optional[bool] = None
    ionic_conv: Optional[bool] = None

