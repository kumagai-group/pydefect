# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from dataclasses import dataclass
from typing import List, Optional

from monty.json import MSONable
from pydefect.input_maker.append_interstitial import append_interstitial
from pydefect.util.coords import pretty_coords
from pydefect.util.structure_tools import Coordination
from pymatgen.core import Structure
from tabulate import tabulate
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
    threshold_frac: Optional[float]
    threshold_abs: Optional[float]
    min_dist: float
    tol: float
    radius: float


@dataclass
class VolumetricDataLocalExtrema(MSONable, ToJsonFileMixIn):
    unit_cell: Structure
    is_min: bool
    extrema_points: List[CoordInfo]
    info: str
    params: VolumetricDataAnalyzeParams

    def __str__(self):
        min_or_max = "min" if self.is_min else "max"
        extrema_points = [["#", "site_sym", "coordination", "frac_coords", "quantity"]]
        for idx, ep in enumerate(self.extrema_points, 1):
            extrema_points.append([idx,
                                   ep.site_symmetry,
                                   ep.coordination.distance_dict,
                                   pretty_coords(ep.frac_coords[0]),
                                   f"{ep.quantities[0]:.2g}"])
        lines = [f"info: {self.info}",
                 f"min_or_max: {min_or_max}",
                 f"extrema_points:",
                 tabulate(extrema_points, tablefmt="plain")]
        return "\n".join(lines)

    def append_sites_to_supercell_info(self, supercell_info, indices):
        frac_coords, infos = [], []
        for idx, ep in enumerate(self.extrema_points, 1):
            if idx in indices:
                frac_coords.append(ep.frac_coords[0])
                infos.append(f"{self.info} #{idx}")

        result = append_interstitial(supercell_info,
                                     self.unit_cell,
                                     frac_coords=frac_coords,
                                     infos=infos)
        return result

