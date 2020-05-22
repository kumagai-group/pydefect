# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass

import numpy as np
from monty.json import MSONable


@dataclass
class Unitcell(MSONable):
    vbm: float
    cbm: float
    ele_dielectric_const: np.array
    ion_dielectric_const: np.array

    @property
    def dielectric_constant(self):
        return self.ele_dielectric_const + self.ion_dielectric_const


