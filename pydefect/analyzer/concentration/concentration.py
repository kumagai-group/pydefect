# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from dataclasses import dataclass
from typing import List, Optional, Dict

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
        return self.concentrations

    def __str__(self):
        return tabulate([self.str_list], headers=self.str_header, floatfmt=".2e")

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
    """Concentration per cell at a given Fermi level."""
    Ef: float
    carrier: CarrierConcentration
    defects: List[DefectConcentration]

    @property
    def pinning_levels(self) -> List[float]:
        result = [float("-inf"), float("inf")]
        for d in self.defects:
            lower, upper = d.pinning_levels
            result = [max([result[0], lower]), min([result[1], upper])]
        return result

    @property
    def str_header(self):
        result = ["Ef", "p", "n"]
        for d in self.defects:
            result.extend(d.str_header)
        result.extend(["net charge", "net ratio"])
        return result

    @property
    def str_list(self):
        result = [self.Ef, self.carrier.p, self.carrier.n]
        for d in self.defects:
            result.extend(d.str_list)
        result.extend([self.net_charge, self.net_abs_ratio])
        return result

    def __str__(self):
        str_list = [self.str_list]
        floatfmt = [".3f"] + [".1e"] * (len(str_list[0]) - 1)
        return tabulate(str_list, headers=self.str_header, floatfmt=floatfmt)

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
    # Ex: pinning_levels[specie_name] = [float("-inf"), 1.0]
    # The VBM is set to zero.
    pinning_levels: Dict[str, List[Optional[float]]] = None
    equilibrium_concentration: Optional[Concentration] = None
    T_before_quench: float = None

    def __post_init__(self):
        self.concentrations.sort(key=lambda x: x.Ef)

    @property
    def pinning_level(self):
        if self.pinning_levels:
            lower = max(i[0] for i in self.pinning_levels.values())
            upper = min(i[1] for i in self.pinning_levels.values())
            return [lower, upper]

    def __str__(self):
        result = [f"T: {self.T}"]

        table = []
        if self.pinning_levels:
            for k, v in self.pinning_levels.items():
                table.append([k, v[0], v[1]])
            result.append("Pinning levels:")
            result.append(tabulate(table, floatfmt=".3f"))

        table = []
        for c in self.concentrations:
            table.append(c.str_list)
        floatfmt = [".3f"] + [".1e"] * (len(c.str_list) - 1)
        result.append(tabulate(table, c.str_header, floatfmt=floatfmt))
        result.append("-"*50)
        if self.equilibrium_concentration:
            result.append("Equilibrium concentration")
            result.append(self.equilibrium_concentration.__str__())
        return "\n".join(result)

    @property
    def most_neutral_concentration(self) -> Concentration:
        return min(self.concentrations, key=lambda x: abs(x.net_charge))

    def next_Ef_to_neutral_concentration(
            self, n_Ef: int = 10) -> Optional[List[float]]:
        """
        Return multiple concentrations that are closest to a charge neutral
        condition along the negative and positive net charges.

        n_Ef (int): Number of the sampled Fermi level in between end points.

        :return:
            List of the Fermi level (float)
        """
        try:
            c_minus = min([c for c in self.concentrations if c.net_charge > 0],
                          key=lambda x: x.net_charge)
            c_plus = min([c for c in self.concentrations if c.net_charge <= 0],
                         key=lambda x: abs(x.net_charge))
        except ValueError:
            return None

        return np.linspace(
            c_minus.Ef, c_plus.Ef, n_Ef + 2, endpoint=True).tolist()

    # def check_end_point(self, c_minus, c_plus):
    #     end_points = [self.concentrations[0], self.concentrations[-1]]
    #     return c_minus in end_points or c_plus in end_points


