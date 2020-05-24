# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List

import numpy as np
from monty.json import MSONable

from pydefect.corrections.abstract_correction import Correction


@dataclass
class ExtendedFnvCorrection(Correction, MSONable):
    """
    species: Species except for the defect. e.g., ["Mg", "Mg", ..., "O", ..]
    atomic_coords: Fractional coordinates except for the defect.
    pc_pot (list of float):
        List of point-charge potential from the defect for all the
        atomic sites.
    defect_region_radius (float):
        Maximum radius of a sphere touching to the lattice plane, used
        for defining the outside region of the defect.
    """
    charge: int
    point_charge_correction: float
    defect_region_radius: float
    species: List[str]
    distances: List[float]
    site_ele_pot: List[float]
    point_charge_pot: List[float]
    additional_correction: float = 0.0

    @property
    def diff_pot(self):
        return [e - p for e, p in zip(self.site_ele_pot, self.point_charge_pot)]

    @property
    def ave_pot_diff(self):
        return np.mean([p for p, d in zip(self.diff_pot, self.distances)
                        if d > self.defect_region_radius])

    @property
    def alignment_correction(self) -> float:
        return - self.ave_pot_diff * self.charge

    @property
    def correction_energy(self) -> float:
        return (self.point_charge_correction + self.alignment_correction
                + self.additional_correction)
