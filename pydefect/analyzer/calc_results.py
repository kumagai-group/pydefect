# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Optional

from monty.json import MSONable
from pymatgen import IStructure

from pydefect.mix_in import ToJsonFileMixIn


@dataclass
class CalcResults(MSONable, ToJsonFileMixIn):
    structure: IStructure
    site_symmetry: str
    energy: float
    magnetization: float
    kpoint_coords: List[List[float]]
    kpoint_weights: List[float]
    potentials: List[float]
    vbm_info: float
    cbm_info: float
    fermi_level: float
    electronic_conv: Optional[bool] = None
    ionic_conv: Optional[bool] = None



