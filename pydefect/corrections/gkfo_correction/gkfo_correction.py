# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List

import numpy as np
from pydefect.corrections.abstract_correction import Correction
from pydefect.corrections.efnv_correction.efnv_correction import \
    ExtendedFnvCorrection, PotentialSite


@dataclass
class GkfoCorrection(Correction):
    """
    potential of gkfo_sites is the difference of potential in the final state
    from that in the initial state.
    And, pc_potential is caused by the additional_charge.
    """
    efnv_correction: ExtendedFnvCorrection
    additional_charge: int
    pc_2nd_term: float
    gkfo_sites: List[PotentialSite]
    ave_dielectric_tensor: float
    ave_electronic_dielectric_tensor: float

    @property
    def charge(self):
        return self.efnv_correction.charge

    @property
    def defect_region_radius(self):
        return self.efnv_correction.defect_region_radius

    @property
    def pc_1st_term(self):
        return (self.efnv_correction.point_charge_correction * 2 / self.charge
                * self.additional_charge)

    @property
    def sum_pc_correction(self):
        return self.pc_1st_term + self.pc_2nd_term

    @property
    def average_potential_diff_by_addition(self):
        return np.mean([s.diff_pot for s in self.gkfo_sites
                        if s.distance > self.defect_region_radius])

    @property
    def alignment_1st_term(self) -> float:
        return - self.average_potential_diff_by_addition * self.additional_charge

    @property
    def alignment_2nd_term(self) -> float:
        return - self.efnv_correction.average_potential_diff * self.additional_charge

    @property
    def alignment_3rd_term(self) -> float:
        return - (self.ave_electronic_dielectric_tensor
                  / self.ave_dielectric_tensor
                  * self.charge * self.average_potential_diff_by_addition)

    @property
    def sum_alignment_term(self) -> float:
        return self.alignment_1st_term + self.alignment_2nd_term + self.alignment_3rd_term

    @property
    def correction_energy(self) -> float:
        return self.sum_pc_correction + self.sum_alignment_term


