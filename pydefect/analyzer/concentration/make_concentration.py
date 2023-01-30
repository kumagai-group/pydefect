# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from abc import ABCMeta, abstractmethod
from dataclasses import dataclass
from typing import Dict, List

import numpy as np
from pydefect.analyzer.concentration.concentration import DefectConcentration, CarrierConcentration, Concentration, \
    ConcentrationByFermiLevel
from pydefect.analyzer.concentration.distribution_function import \
    boltzmann_dist, fermi_dirac
from pydefect.analyzer.defect_energy import ChargeEnergies, SingleChargeEnergies
from vise.analyzer.dos_data import DosData


class MakeCarrierConcentrations:
    def __init__(self,
                 dos_data: DosData,
                 charge_energies: ChargeEnergies,
                 degeneracy: Dict[str, Dict[int, int]],
                 T: float):
#                 V: float):
        self._fermi_level = np.mean(dos_data.vertical_lines)
        self._T = T
        self.charge_energies = charge_energies.charge_energies_dict
        self._degeneracy = degeneracy
        self._tdos = dos_data.total[0]
        self._make_vb_dos(dos_data)
        self._make_cb_dos(dos_data)

    def _make_vb_dos(self, dos_data):
        energy_range = dos_data.energies < self._fermi_level
        energies = np.array(dos_data.energies)[energy_range].tolist()
        doses = self._tdos[energy_range].tolist()
        self.vb_dos = VBDos(energies, doses)

    def _make_cb_dos(self, dos_data):
        energy_range = dos_data.energies > self._fermi_level
        energies = np.array(dos_data.energies)[energy_range].tolist()
        doses = self._tdos[energy_range].tolist()
        self.cb_dos = CBDos(energies, doses)

    def _make_defect_concentration(self,
                                   name: str,
                                   single_energies: SingleChargeEnergies):
        charges, concentrations = [], []
        for (c, e) in single_energies.charge_energies:
            charges.append(c)
            deg = self._degeneracy[name][c]
            concentrations.append(boltzmann_dist(e, self._T) * deg)

        return DefectConcentration(name, charges, concentrations)

    def _make_all_concentration(self, Ef):
        p = self.vb_dos.carrier_concentration(Ef, self._T)
        n = self.cb_dos.carrier_concentration(Ef, self._T)
        carrier = CarrierConcentration(p, n)

        defects = []
        for name, single in self.charge_energies.items():
            concentration = self._make_defect_concentration(name, single)
            defects.append(concentration)

        return Concentration(Ef, carrier, defects)

    def make_concentrations_by_fermi_level(self, Efs: List[float]
                                           ) -> ConcentrationByFermiLevel:
        concentrations = [self._make_all_concentration(Ef) for Ef in Efs]
        return ConcentrationByFermiLevel(self._T, concentrations)


def equilibrium_concentration(make_cc: MakeCarrierConcentrations,
                              e_min: float,
                              e_max: float,
                              n: int = 10,
                              n_iter: int = 10,
                              net_abs_ratio: float = 1.0e-5) -> Concentration:
    Efs = np.linspace(e_min, e_max, n + 1).tolist()

    for _ in range(n_iter):
        cons_Ef = make_cc.make_concentrations_by_fermi_level(Efs)
        if cons_Ef.most_neutral_concentration.net_abs_ratio < net_abs_ratio:
            return cons_Ef.most_neutral_concentration
        Efs = cons_Ef.next_Ef_to_neutral_concentration(n)

    raise ValueError("No convergence is obtained.")


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