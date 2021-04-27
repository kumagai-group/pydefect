# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from dataclasses import dataclass
from typing import List

from monty.json import MSONable
from pydefect.util.structure_tools import Coordination
from pymatgen import Structure
from vise.util.mix_in import ToJsonFileMixIn
from vise.util.typing import Coords


@dataclass
class CoordInfo(MSONable):
    frac_coords: Coords
    site_symmetry: str
    coordination: Coordination
    quantity: float


@dataclass
class LocalExtrema(MSONable, ToJsonFileMixIn):
    unit_cell: Structure
    is_min: bool
    extrema_points: List[CoordInfo]
    info: str
