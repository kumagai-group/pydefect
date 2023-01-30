# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.

import pytest
from numpy.testing import assert_almost_equal
from pydefect.analyzer.concentration.concentration import CarrierConcentration, \
    DefectConcentration, Concentration, ConcentrationByFermiLevel


def test_carrier_concentration():
    concentration = CarrierConcentration(p=1.23*10**15, n=0.23*10**15)
    assert_almost_equal(concentration.net_charge, 1.0*10**15)


def test_defect_concentration():
    concentration = DefectConcentration(name="Va_O1",
                                        charges=[0, 1, 2],
                                        concentrations=[1.0, 2.0, 3.0])
    assert_almost_equal(concentration.net_charge, 8.0)
    assert_almost_equal(concentration.total_concentration, 6.0)


@pytest.fixture
def concentration():
    return Concentration(
        Ef=0.0,
        carrier=CarrierConcentration(p=1.23, n=0.23),
        defects=[DefectConcentration(name="Va_O1",
                                     charges=[0, 1, 2],
                                     concentrations=[1.0, 2.0, 3.0])])


def test_concentrations(concentration):
    assert concentration.net_charge == 9.0
    assert concentration.abs_charge == 9.46
    assert concentration.net_abs_ratio == 9.0 / 9.46


def test_concentration_by_fermi_level(mocker):
    mock_1 = mocker.Mock(spec=Concentration, autospec=True)
    mock_2 = mocker.Mock(spec=Concentration, autospec=True)
    mock_3 = mocker.Mock(spec=Concentration, autospec=True)
    mock_4 = mocker.Mock(spec=Concentration, autospec=True)

    mock_1.net_charge = 1.0
    mock_2.net_charge = 0.1
    mock_3.net_charge = -1.0
    mock_4.net_charge = -2.0

    mock_1.Ef = 0.0
    mock_2.Ef = 1.0
    mock_3.Ef = 2.0
    mock_4.Ef = 3.0

    con = ConcentrationByFermiLevel(
        T=300, concentrations=[mock_1, mock_2, mock_3, mock_4])
    assert con.most_neutral_concentration == mock_2
    assert con.next_Ef_to_neutral_concentration(3) == [1.25, 1.5, 1.75]


