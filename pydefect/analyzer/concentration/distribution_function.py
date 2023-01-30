# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from numpy import exp
from scipy.constants import physical_constants


k = physical_constants["Boltzmann constant in eV/K"][0]


def fermi_dirac(delta_E: float, T: float):  # delta_E is in eV.
    return 1. / (exp(delta_E / (k * T)) + 1.)


def boltzmann_dist(delta_E: float, T: float):
    return exp(- delta_E / (k * T))

