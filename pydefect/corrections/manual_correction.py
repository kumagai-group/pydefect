# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pydefect.corrections.abstract_correction import Correction


class NoCorrection(Correction):
    @property
    def correction_energy(self):
        return 0.0


class ManualCorrection(Correction):
    def __init__(self, correction_energy: float):
        self._correction_energy = correction_energy

    @property
    def correction_energy(self):
        return self._correction_energy
