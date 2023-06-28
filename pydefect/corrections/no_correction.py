# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from pydefect.corrections.abstract_correction import Correction


class NoCorrection(Correction):
    @property
    def correction_energy(self) -> float:
        return 0.0

    @property
    def correction_dict(self) -> dict:
        return {}
