# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Optional

from monty.json import MSONable
from pymatgen import IStructure
from vise.analyzer.band_edge_properties import BandEdge
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class CalcResults(MSONable, ToJsonFileMixIn):
    structure: IStructure
    site_symmetry: str
    energy: float
    magnetization: float
    kpoint_coords: List[List[float]]
    kpoint_weights: List[float]
    # potential acting on the positive unit charge, whose sign is reserved from
    # vasp convention of atomic site potential.
    potentials: List[float]
    vbm_info: BandEdge
    cbm_info: BandEdge
    fermi_level: float
    electronic_conv: Optional[bool] = None
    ionic_conv: Optional[bool] = None

    @property
    def vbm(self):
        return self.vbm_info.energy

    @property
    def cbm(self):
        return self.cbm_info.energy


