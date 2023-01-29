# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from abc import ABCMeta, abstractmethod
from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
from matplotlib.axes import Axes
from numpy import exp
from pydefect.analyzer.defect_energy import SingleChargeEnergies, ChargeEnergies
from scipy.constants import  physical_constants
from vise.analyzer.dos_data import DosData

k = physical_constants["Boltzmann constant in eV/K"][0]


@dataclass
class CarrierConcentration:
    """Concentration per cell"""
    p: float
    n: float

    @property
    def net_charge(self):
        return self.p - self.n


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


@dataclass
class ConcentrationByFermiLevel:
    """Concentration per cell"""
    T: float
    concentrations: List[Concentration]
    quench_from_T: float = None

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


def fermi_dirac(delta_E: float, T: float):  # delta_E is in eV.
    return 1. / (exp(delta_E / (k * T)) + 1.)


def boltzmann_dist(delta_E: float, T: float):
    return exp(- delta_E / (k * T))


class MakeCarrierConcentrations:
    def __init__(self,
                 dos_data: DosData,
                 charge_energies: ChargeEnergies,
                 Efs: List[float],
                 T: float):
        self._fermi_level = np.mean(dos_data.vertical_lines)
        self.T = T
        self.charge_energies = charge_energies.charge_energies_dict
        self.e_min = charge_energies.e_min
        self.e_max = charge_energies.e_max

        self._tdos = dos_data.total[0]
        self._make_vb_dos(dos_data)
        self._make_cb_dos(dos_data)

        self.cons_by_Ef = self._make_concentrations_by_fermi_level(Efs)

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
            concentrations.append(boltzmann_dist(e, self.T))
        return DefectConcentration(name, charges, concentrations)

    def _make_all_concentration(self, Ef):
        p = self.vb_dos.carrier_concentration(Ef, self.T)
        n = self.cb_dos.carrier_concentration(Ef, self.T)
        carrier = CarrierConcentration(p, n)

        defects = []
        for name, single in self.charge_energies.items():
            concentration = self._make_defect_concentration(name, single)
            defects.append(concentration)

        return Concentration(Ef, carrier, defects)

    def _make_concentrations_by_fermi_level(self, Efs: List[float]):
        concentrations = [self._make_all_concentration(Ef) for Ef in Efs]
        return ConcentrationByFermiLevel(self.T, concentrations)


def plot_pn(cc: ConcentrationByFermiLevel, ax: Axes):
    Efs, ps, ns = [], [], []
    for cc in cc.concentrations:
        Efs.append(cc.Ef)
        ps.append(cc.carrier.p)
        ns.append(cc.carrier.n)

    ax.set_yscale("log")
    ax.plot(Efs, ps)
    ax.plot(Efs, ns)


def plot_defect_concentration(cc: ConcentrationByFermiLevel, ax: Axes):
    Efs, ddd = [], defaultdict(list)
    for cc in cc.concentrations:
        Efs.append(cc.Ef)
        for defect in cc.defects:
            ddd[defect.name].append(defect.total_concentration)

    ax.set_yscale("log")
    for dd in ddd:
        ax.plot(Efs, dd)


# # change to class
# def make_carrier_concentrations(dos_data: DosData, T: float):
#     dos_interval = dos_data.energies[1] - dos_data.energies[0]
#     tdos = dos_data.total[0]
#     fermi_level = np.mean(dos_data.vertical_lines)
#     vb_dos = tdos[dos_data.energies < fermi_level]
#     cb_dos = tdos[dos_data.energies > fermi_level]

    # return