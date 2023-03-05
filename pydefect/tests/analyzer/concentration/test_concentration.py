# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.

import pytest
from numpy.testing import assert_almost_equal
from pydefect.analyzer.concentration.concentration import \
    CarrierConcentration, DefectConcentration, Concentration, \
    ConcentrationByFermiLevel


@pytest.fixture
def carrier_concentration():
    return CarrierConcentration(p=1.23*10**15, n=0.23*10**15)


@pytest.fixture
def defect_concentration():
    return DefectConcentration(name="Va_Mg1",
                               charges=[-1, 0, 1, 2],
                               concentrations=[1.1, 1.0, 2.0, 3.0])


@pytest.fixture
def concentration():
    return Concentration(
        Ef=0.0,
        carrier=CarrierConcentration(p=1.23, n=0.23),
        defects=[DefectConcentration(name="Va_O1",
                                     charges=[0, 1, 2],
                                     concentrations=[1.0, 2.0, 3.0]),
                 DefectConcentration(name="Va_Mg1",
                                     charges=[0, 1, 2],
                                     concentrations=[1.0, 2.0, 3.0]),
                 ])


def test_carrier_concentration_net_charge(carrier_concentration):
    assert_almost_equal(carrier_concentration.net_charge, 1.0*10**15)


def test_carrier_concentration_str(carrier_concentration):
    print(carrier_concentration)


def test_defect_concentration(defect_concentration):
    assert_almost_equal(defect_concentration.net_charge, 6.9)
    assert_almost_equal(defect_concentration.abs_charge, 9.1)
    assert_almost_equal(defect_concentration.total_concentration, 7.1)


def test_defect_concentration_str(defect_concentration):
    print(defect_concentration)


def test_concentrations(concentration):
    assert concentration.net_charge == 17.0
    assert concentration.abs_charge == 17.46
    assert concentration.net_abs_ratio == 17.0 / 17.46


def test_concentration_str(concentration):
    print(concentration)


@pytest.fixture
def con_by_Ef(concentration):
    va_o1 = DefectConcentration("Va_O1", [0], [1.0])
    c1 = Concentration(0.0, CarrierConcentration(1.0, 0.0), [va_o1])
    c2 = Concentration(1.0, CarrierConcentration(0.1, 0.0), [va_o1])
    c3 = Concentration(2.0, CarrierConcentration(0.0, 1.0), [va_o1])
    c4 = Concentration(3.0, CarrierConcentration(0.0, 2.0), [va_o1])
    return ConcentrationByFermiLevel(
        T=300, concentrations=[c1, c2, c3, c4],
        equilibrium_concentration=concentration,
        pinning_levels={"Va_O1": [float("-inf"), 1.0]})


def test_concentration_by_fermi_level(con_by_Ef):
    assert con_by_Ef.most_neutral_concentration == con_by_Ef.concentrations[1]


def test_next_Ef_to_neutral_concentration(con_by_Ef):
    actual = con_by_Ef.next_Ef_to_neutral_concentration(3)
    assert actual == [1.0, 1.25, 1.5, 1.75, 2.0]


def test_pinning_level(concentration):
    con_by_Ef = ConcentrationByFermiLevel(
        T=1, concentrations=[],
        equilibrium_concentration=concentration,
        pinning_levels={"Va_O1": [float("-inf"), 1.0],
                        "Va_Mg1": [-1.0, 0.5]})
    assert con_by_Ef.pinning_level == [-1.0, 0.5]


def test_con_by_Ef_str(con_by_Ef):
    print(con_by_Ef)
