# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from dataclasses import dataclass
from typing import List

import numpy as np


@dataclass
class CarrierConcentration:
    """Concentration per cell"""
    p: float
    n: float

    @property
    def net_charge(self):
        return self.p - self.n

    @property
    def abs_charge(self):
        return self.p + self.n


@dataclass
class DefectConcentration:
    """Concentration per cell"""
    name: str
    charges: List[int]
    concentrations: List[float]

    @property
    def net_charge(self):
        return sum(c * con for c, con in zip(self.charges, self.concentrations))

    @property
    def abs_charge(self):
        return sum(abs(c * con)
                   for c, con in zip(self.charges, self.concentrations))

    @property
    def total_concentration(self):
        return sum(self.concentrations)


@dataclass
class Concentration:
    """Concentration per cell"""
    Ef: float
    carrier: CarrierConcentration
    defects: List[DefectConcentration]

    @property
    def net_charge(self):
        total_defect_charge = sum(d.net_charge for d in self.defects)
        return self.carrier.net_charge + total_defect_charge

    @property
    def abs_charge(self):
        defect_abs_charge = sum(d.abs_charge for d in self.defects)
        return self.carrier.abs_charge + defect_abs_charge

    @property
    def net_abs_ratio(self):
        return abs(self.net_charge) / self.abs_charge


@dataclass
class ConcentrationByFermiLevel:
    """Concentration per cell"""
    T: float
    concentrations: List[Concentration]
#    quench_from_T: float = None

    def __post_init__(self):
        self.concentrations.sort(key=lambda x: x.Ef)

    @property
    def most_neutral_concentration(self) -> Concentration:
        return min(self.concentrations, key=lambda x: abs(x.net_charge))

    def next_Ef_to_neutral_concentration(self, n_Ef: int = 10) -> List[float]:
        """
        Return two concentrations that are closest to a charge neutral condition
        along the negative and positive net charges.

        :return:
        """
        c_minus = min([c for c in self.concentrations if c.net_charge > 0],
                      key=lambda x: x.net_charge)
        c_plus = min([c for c in self.concentrations if c.net_charge <= 0],
                     key=lambda x: abs(x.net_charge))

        return np.linspace(
            c_minus.Ef, c_plus.Ef, n_Ef + 1, endpoint=False)[1:].tolist()


