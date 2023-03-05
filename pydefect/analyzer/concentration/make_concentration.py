# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from abc import ABCMeta, abstractmethod
from dataclasses import dataclass
from typing import List, Dict, Optional

import numpy as np
from monty.json import MSONable
from pydefect.analyzer.concentration.concentration import DefectConcentration, \
    CarrierConcentration, Concentration, \
    ConcentrationByFermiLevel
from pydefect.analyzer.concentration.degeneracy import Degeneracies
from pydefect.analyzer.concentration.distribution_function import \
    boltzmann_dist, fermi_dirac
from pydefect.analyzer.defect_energy import ChargeEnergies, SingleChargeEnergies
from tabulate import tabulate
from vise.util.logger import get_logger
from vise.util.mix_in import ToJsonFileMixIn


logger = get_logger(__name__)


@dataclass
class TotalDos(MSONable, ToJsonFileMixIn):
    # the energies are in absolute scale.
    energies: List[float]  # VBM is set to 0.
    dos: List[float]
    volume: float  # in Å^3
    vbm: float
    cbm: float

    @property
    def fermi_level(self):
        return (self.vbm + self.cbm) / 2

    def __str__(self):
        result = [f"Fermi level: {self.fermi_level}",
                  f"VBM: {self.vbm}",
                  f"CBM: {self.cbm}",
                  f"Volume: {self.volume}"]
        list_ = []
        for e, d in zip(self.energies, self.dos):
            list_.append([e, d])
        result.append(tabulate(list_, headers=["energy", "dos"],
                               floatfmt=".2e"))
        return "\n".join(result)


class MakeCarrierConcentrations:
    def __init__(self, total_dos: TotalDos, T: float):
        self._fermi_level = total_dos.fermi_level
        self._tdos = np.array(total_dos.dos)
        self._abs_energies = np.array(total_dos.energies)
        self._V = total_dos.volume / 10**24  # in cm^3
        self._vbm = total_dos.vbm
        self.T = T
        self._make_vb_dos()
        self._make_cb_dos()

    def _make_vb_dos(self):
        energy_range = self._abs_energies < self._fermi_level
        energies = np.array(self._abs_energies)[energy_range]
        doses = self._tdos[energy_range].tolist()
        self.vb_dos = VBDos((energies - self._vbm).tolist(), doses)

    def _make_cb_dos(self):
        energy_range = self._abs_energies > self._fermi_level
        energies = np.array(self._abs_energies)[energy_range]
        doses = self._tdos[energy_range].tolist()
        self.cb_dos = CBDos((energies - self._vbm).tolist(), doses)

    def carrier_concentration(self, Ef):
        p = self.vb_dos.carrier_concentration(Ef, self.T) / self._V
        n = self.cb_dos.carrier_concentration(Ef, self.T) / self._V
        return CarrierConcentration(p, n)

    def _concentration(self, Ef):
        return Concentration(Ef, self.carrier_concentration(Ef), [])

    def make_concentrations_by_fermi_level(self, Efs: List[float],
                                           ) -> ConcentrationByFermiLevel:
        concentrations = [self._concentration(Ef) for Ef in Efs]
        return ConcentrationByFermiLevel(self.T, concentrations, None)


class MakeConcentrations:
    def __init__(self,
                 total_dos: TotalDos,
                 charge_energies: ChargeEnergies,
                 degeneracies: Degeneracies,
                 T: float,
                 fixed_defect_concentrations: Dict[str, float] = None):

        self._tdos = np.array(total_dos.dos)
        self._abs_energies = np.array(total_dos.energies)
        self._V = total_dos.volume / 10**24  # in cm^3
        self.fixed_con = fixed_defect_concentrations

        self.make_carrier_concentrations = MakeCarrierConcentrations(total_dos, T)

        self.T = T
        self.charge_energies = charge_energies
        self._degeneracies = degeneracies

    def _make_defect_concentration(self,
                                   name: str,
                                   Ef: float,
                                   single_energies: SingleChargeEnergies):
        charges, concentrations = [], []
        for (c, e) in single_energies.charge_energies_at_ef(Ef):
            charges.append(c)
            deg = self._degeneracies[name][c].degeneracy
            concentrations.append(boltzmann_dist(e, self.T) * deg / self._V)

        concentrations = self._redistribute(concentrations, name)
        return DefectConcentration(name, charges, concentrations)

    def _calc_pinning(self, single_energies):
        pin_level = single_energies.pinning_level(float("-inf"), float("inf"))
        try:
            lower = pin_level[0][0]
        except TypeError:
            lower = float("-inf")
        try:
            upper = pin_level[1][0]
        except TypeError:
            upper = float("inf")
        return [lower, upper]

    def _redistribute(self, concentrations, name):
        if self.fixed_con:
            return redistribute_concentration(
                concentrations, self.fixed_con[name])
        return concentrations

    def _make_all_concentration(self, Ef: float):
        carrier = self.make_carrier_concentrations.carrier_concentration(Ef)

        defects = []
        for name, single in self.charge_energies.charge_energies_dict.items():
            concentration = self._make_defect_concentration(name, Ef, single)
            defects.append(concentration)

        return Concentration(Ef, carrier, defects)

    def make_concentrations_by_fermi_level(self, Efs: List[float],
                                           ) -> ConcentrationByFermiLevel:
        concentrations = [self._make_all_concentration(Ef) for Ef in Efs]
        pinning_levels = self.charge_energies.pinning_levels
        return ConcentrationByFermiLevel(self.T,
                                         concentrations,
                                         pinning_levels)




def redistribute_concentration(
        concentrations: List[float], total: float) -> List[float]:
    factor = total / np.sum(concentrations)
    return (np.array(concentrations) * factor).tolist()


def equilibrium_concentration(make_cc: MakeConcentrations,
                              e_min: float,
                              e_max: float,
                              n: int = 10,
                              n_iter: int = 10,
                              net_abs_ratio: float = 1.0e-5
                              ) -> Optional[Concentration]:
    Efs = np.linspace(e_min, e_max, n + 1).tolist()

    for i in range(n_iter):
        logger.info(f"Calc equilibrium concentration: iteration {i}")
        cons_Ef = make_cc.make_concentrations_by_fermi_level(Efs)

        if cons_Ef.most_neutral_concentration.net_abs_ratio < net_abs_ratio:
            logger.info(f"Equilibrium concentration is found.")
            return cons_Ef.most_neutral_concentration

        Efs = cons_Ef.next_Ef_to_neutral_concentration(n)

        if Efs is None:
            logger.warning(f"Charge balance is out of {e_min}--{e_max}.")
            return

    logger.warning("No convergence is obtained.")


@dataclass
class Dos(metaclass=ABCMeta):
    energies: List[float]
    doses: List[float]

    def carrier_concentration(self, Ef, T) -> float:
        result = 0.0
        for E, dos in zip(self.energies, self.doses):
            result += self.interval * dos * self._fermi_dirac(Ef, E, T)
        return result

    @staticmethod
    @abstractmethod
    def _fermi_dirac(Ef, E, T):
        pass

    @property
    def interval(self) -> float:
        return self.energies[1] - self.energies[0]


class VBDos(Dos):
    carrier_type = "p"

    @staticmethod
    def _fermi_dirac(Ef, E, T):
        return fermi_dirac(Ef - E, T)


class CBDos(Dos):
    carrier_type = "n"

    @staticmethod
    def _fermi_dirac(Ef, E, T):
        return fermi_dirac(E - Ef, T)