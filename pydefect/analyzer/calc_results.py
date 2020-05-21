# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from monty.json import MSONable
from monty.serialization import loadfn
from pymatgen import IStructure


@dataclass
class CalcResults(MSONable):
    structure: IStructure
    site_symmetry: str
    energy: float
    magnetization: float
    kpoint_coords: List[List[float]]
    kpoint_weights: List[float]
    electrostatic_potential: List[float]
    vbm_info: float
    cbm_info: float
    fermi_level: float
    electronic_conv: Optional[bool] = None
    ionic_conv: Optional[bool] = None

    def to_json_file(self, filename: str = "calc_results.json") -> None:
        Path(filename).write_text(self.to_json())

    @classmethod
    def from_json_file(cls, filename: str = "calc_results.json"):
        return loadfn(filename)
