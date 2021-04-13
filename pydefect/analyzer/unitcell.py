# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List

import numpy as np
from monty.json import MSONable
from monty.serialization import loadfn
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

    def to_yaml(self, filename: str = "unitcell.yaml") -> None:
        text = f"""vbm: {self.vbm}
cbm: {self.cbm}
ele_dielectric_const:
- {self.ele_dielectric_const[0]}
- {self.ele_dielectric_const[1]}
- {self.ele_dielectric_const[2]}
ion_dielectric_const:
- {self.ion_dielectric_const[0]}
- {self.ion_dielectric_const[1]}
- {self.ion_dielectric_const[2]}
"""
        with open(filename, 'w') as file:
            file.write(text)

    @classmethod
    def from_yaml(cls, filename: str = "unitcell.yaml"):
        return cls(**loadfn(filename))
