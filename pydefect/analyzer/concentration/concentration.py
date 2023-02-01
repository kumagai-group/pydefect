# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from monty.json import MSONable
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class CarrierConcentration(MSONable):
    """Concentration per cell"""
    p: float
    n: float

    def __str__(self):
        return f"p: {self.p:.2e}, n: {self.n:.2e}"

    @property
    def net_charge(self):
        return self.p - self.n

    @property
    def abs_charge(self):
        return self.p + self.n


@dataclass
class DefectConcentration(MSONable):
    """Concentration per cell"""
    name: str
    charges: List[int]
    concentrations: List[float]

    @property
    def str_header(self):
        return [f"{self.name}_{c}" for c in self.charges]

    @property
    def str_list(self):
        return [f"{con:.2e}" for con in self.concentrations]

    def __str__(self):
        return tabulate([self.str_list], headers=self.str_header)

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
class Concentration(MSONable):
    """Concentration per cell"""
    Ef: float
    carrier: CarrierConcentration
    defects: List[DefectConcentration]

    @property
    def str_header(self):
        result = ["Ef", "p", "n"]
        for d in self.defects:
            result.extend(d.str_header)
        result.extend(["net charge", "net ratio"])
        return result

    @property
    def str_list(self):
        result = [f"{self.Ef:.2f}",
                  f"{self.carrier.p:.2e}", f"{self.carrier.n:.2e}"]
        for d in self.defects:
            result.extend(d.str_list)
        result.extend([f"{self.net_charge:.1e}", f"{self.net_abs_ratio:.1e}"])
        return result

    def __str__(self):
        return tabulate([self.str_list], headers=self.str_header)

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
class ConcentrationByFermiLevel(MSONable, ToJsonFileMixIn):
    """Concentration per cell"""
    T: float
    concentrations: List[Concentration]
    equilibrium_concentration: Optional[Concentration]
    T_before_quench: float = None

    def __post_init__(self):
        self.concentrations.sort(key=lambda x: x.Ef)

    def __str__(self):
        result = [f"T: {self.T}"]
        table = []
        for c in self.concentrations:
            table.append(c.str_list)
        result.append(tabulate(table, c.str_header))
        result.append("-"*50)
        if self.equilibrium_concentration:
            result.append("Equilibrium concentration")
            result.append(self.equilibrium_concentration.__str__())
        return "\n".join(result)

    @property
    def most_neutral_concentration(self) -> Concentration:
        return min(self.concentrations, key=lambda x: abs(x.net_charge))

    def next_Ef_to_neutral_concentration(self, n_Ef: int = 10
                                         ) -> Optional[List[float]]:
        """
        Return multiple concentrations that are closest to a charge neutral
        condition along the negative and positive net charges.

        :return:
        """
        try:
            c_minus = min([c for c in self.concentrations if c.net_charge > 0],
                          key=lambda x: x.net_charge)
            c_plus = min([c for c in self.concentrations if c.net_charge <= 0],
                         key=lambda x: abs(x.net_charge))
        except ValueError:
            return None

        if self.check_end_point(c_minus, c_plus):
            return None

        return np.linspace(
            c_minus.Ef, c_plus.Ef, n_Ef + 1, endpoint=False)[1:].tolist()

    def check_end_point(self, c_minus, c_plus):
        end_points = [self.concentrations[0], self.concentrations[-1]]
        return c_minus in end_points or c_plus in end_points


