# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from math import exp

import numpy as np
import pytest
from matplotlib import pyplot as plt
from numpy.testing import assert_almost_equal
from pydefect.analyzer.concentration import CarrierConcentration, \
    DefectConcentration, Concentration, ConcentrationByFermiLevel, \
    VBDos, fermi_dirac, CBDos, MakeCarrierConcentrations, boltzmann_dist, \
    plot_pn, plot_defect_concentration
from pydefect.analyzer.defect_energy import SingleChargeEnergies, ChargeEnergies
from vise.analyzer.dos_data import DosData


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
def concentrations():
    return Concentration(
        Ef=0.0,
        carrier=CarrierConcentration(p=1.23, n=0.23),
        defects=[DefectConcentration(name="Va_O1",
                                     charges=[0, 1, 2],
                                     concentrations=[1.0, 2.0, 3.0])])


def test_concentrations(concentrations):
    assert_almost_equal(concentrations.net_charge, 9.0)


def test_concentration_by_fermi_level(mocker, concentrations):
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


def test_fermi_dirac():
    assert_almost_equal(fermi_dirac(0, 1000), 0.5)
    assert_almost_equal(fermi_dirac(1.0, 1 / 8.617333262e-05), 1 / (exp(1) + 1))


def test_boltzmann_dist():
    assert_almost_equal(boltzmann_dist(0, 1000), 1.0)


vb_dos = VBDos(energies=[-1.0, 0.0], doses=[1.0, 1.0])
cb_dos = CBDos(energies=[1.0, 2.0], doses=[1.0, 1.0])


def test_dos_carrier_type():
    assert vb_dos.carrier_type == "p"
    assert cb_dos.carrier_type == "n"


def test_vb_dos():
    actual = vb_dos.carrier_concentration(0.0, 1000)
    expected = 0.5 + fermi_dirac(1.0, 1000)
    assert actual == expected


def test_cb_dos():
    actual = cb_dos.carrier_concentration(1.0, 1000)
    expected = 0.5 + fermi_dirac(1.0, 1000)
    assert actual == expected


single_energies = SingleChargeEnergies([(0, 0.0)])

#
# def test_make_defect_concentration():
#     actual = make_defect_concentration("Va_O1", single_energies, 1000.0)
#     expected = DefectConcentration(name="Va_O1",
#                                    charges=[0],
#                                    concentrations=[1.0])
#     assert actual == expected


def test_make_carrier_concentrations():
    dos_data = DosData(energies=list(np.linspace(-3.0, 3.0, 7)),
                       total=np.array([[1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0]]),
                       pdos=[],
                       vertical_lines=[0.0, 1.5])
    charge_energies = ChargeEnergies({"Va_O1": single_energies},
                                     e_min=0.0,
                                     e_max=1.0)

    make_concentration = MakeCarrierConcentrations(dos_data, charge_energies,
                                                   T=1000.0,
                                                   Efs=[1.0])
    assert make_concentration.vb_dos == VBDos(energies=[-3.0, -2.0, -1.0, 0.0],
                                              doses=[1.0, 1.0, 1.0, 1.0])
    assert make_concentration.cb_dos == CBDos(energies=[1.0, 2.0, 3.0],
                                              doses=[0.0, 1.0, 1.0])

    actual = make_concentration.cons_by_Ef
    expected = ConcentrationByFermiLevel(
        T=1000.0,
        concentrations=[Concentration(Ef=1.0, carrier=CarrierConcentration(p=9.124767652322644e-06, n=9.124767651562896e-06), defects=[DefectConcentration(name='Va_O1', charges=[0], concentrations=[1.0])])])
    assert actual == expected


def test_plot_concentration():
    cc = [Concentration(x, CarrierConcentration(10**(-x), 10**x), [])
          for x in range(1, 10)]
    c = ConcentrationByFermiLevel(300, cc)
    _, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    plot_pn(c, ax1)

    plot_defect_concentration(c, ax2)

    plt.show()


