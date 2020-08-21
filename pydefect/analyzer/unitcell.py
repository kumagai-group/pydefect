# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List

import numpy as np
from monty.json import MSONable
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class Unitcell(MSONable, ToJsonFileMixIn):
    vbm: float
    cbm: float
    ele_dielectric_const: List[List[float]]
    ion_dielectric_const: List[List[float]]

    @property
    def dielectric_constant(self):
        total = (np.array(self.ele_dielectric_const)
                 + np.array(self.ion_dielectric_const))
        return total.tolist()

