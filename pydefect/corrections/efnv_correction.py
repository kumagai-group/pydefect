# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
from monty.json import MSONable
from pydefect.corrections.abstract_correction import Correction
from tabulate import tabulate


@dataclass
class ExtendedFnvCorrection(Correction):
    """
    species: Species except for the defect. e.g., ["Mg", "Mg", ..., "O", ..]
    atomic_coords: Fractional coordinates except for the defect.
    pc_pot (list of float):
        List of point-charge potential from the defect for all the
        atomic sites.
    defect_region_radius (float):
        Maximum radius of a sphere touching to the lattice plane, used
        for defining the outside region of the defect.

    Add units of length and potential
    """
    charge: int
    point_charge_correction: float
    defect_region_radius: float
    sites: List["PotentialSite"]
    defect_coords: Tuple[float, float, float]

    def __str__(self):
        d = [["charge", self.charge],
             ["pc term", self.point_charge_correction],
             ["alignment term", self.alignment_correction],
             ["correction energy", self.correction_energy]]
        return tabulate(d, tablefmt='psql')

    @property
    def average_potential_diff(self):
        return np.mean([s.diff_pot for s in self.sites
                        if s.distance > self.defect_region_radius])

    @property
    def alignment_correction(self) -> float:
        return - self.average_potential_diff * self.charge

    @property
    def correction_energy(self) -> float:
        return self.point_charge_correction + self.alignment_correction

    @property
    def correction_dict(self):
        return {"pc term": self.point_charge_correction,
                "alignment term": self.alignment_correction}


@dataclass
class PotentialSite(MSONable):
    specie: str
    distance: float
    potential: float
    pc_potential: Optional[float]

    @property
    def diff_pot(self):
        return self.potential - self.pc_potential
