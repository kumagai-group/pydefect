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
    site_symmetry: str
    coordination: Coordination
    frac_coords: List[Coords]
    quantities: List[float] = None


@dataclass
class VolumetricDataAnalyzeParams(MSONable):
    threshold_frac: float = None
    threshold_abs: float = None
    min_dist: float = 0.5
    tol: float = 0.5
    radius: float = 0.4


@dataclass
class VolumetricDataLocalExtrema(MSONable, ToJsonFileMixIn):
    unit_cell: Structure
    is_min: bool
    extrema_points: List[CoordInfo]
    info: str
    params: VolumetricDataAnalyzeParams
