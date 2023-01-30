# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from math import exp

from numpy.testing import assert_almost_equal
from pydefect.analyzer.concentration.distribution_function import fermi_dirac, \
    boltzmann_dist


def test_fermi_dirac():
    assert_almost_equal(fermi_dirac(0, 1000), 0.5)
    assert_almost_equal(fermi_dirac(1.0, 1 / 8.617333262e-05), 1 / (exp(1) + 1))


def test_boltzmann_dist():
    assert_almost_equal(boltzmann_dist(0, 1000), 1.0)