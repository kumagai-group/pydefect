# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Dict


@dataclass
class DefectEnergy:
    name: str
    num_atom_change: Dict[str, int]
    charges: List[int]
    energies: List[float]
    corrections: List[float]
    shallow: List[bool]


@dataclass
class DefectEnergies:
    defect_energies: List[DefectEnergy]
    chempot: Dict[str, int]
    keywords: List[str]
    vbm: float
    cbm: float
    supercell_vbm: float
    supercell_cbm: float
