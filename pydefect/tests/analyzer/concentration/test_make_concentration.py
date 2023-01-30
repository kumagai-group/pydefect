# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
import numpy as np
from pydefect.analyzer.concentration.concentration import \
    ConcentrationByFermiLevel, Concentration, CarrierConcentration, \
    DefectConcentration
from pydefect.analyzer.concentration.distribution_function import fermi_dirac
from pydefect.analyzer.concentration.make_concentration import \
    MakeCarrierConcentrations, VBDos, CBDos
from pydefect.analyzer.defect_energy import ChargeEnergies, SingleChargeEnergies
from vise.analyzer.dos_data import DosData


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


def test_make_carrier_concentrations():
    dos_data = DosData(energies=list(np.linspace(-3.0, 3.0, 7)),
                       total=np.array([[1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0]]),
                       pdos=[],
                       vertical_lines=[0.0, 1.5])
    charge_energies = ChargeEnergies({"Va_O1": single_energies},
                                     e_min=0.0,
                                     e_max=1.0)

    degeneracy = {"Va_O1": {0: 1}}

    make_concentration = MakeCarrierConcentrations(
        dos_data=dos_data,
        charge_energies=charge_energies,
        degeneracy=degeneracy,
        T=1000.0)

    assert make_concentration.vb_dos == VBDos(energies=[-3.0, -2.0, -1.0, 0.0],
                                              doses=[1.0, 1.0, 1.0, 1.0])
    assert make_concentration.cb_dos == CBDos(energies=[1.0, 2.0, 3.0],
                                              doses=[0.0, 1.0, 1.0])

    actual = make_concentration.make_concentrations_by_fermi_level(Efs=[1.0])
    expected = ConcentrationByFermiLevel(
        T=1000.0,
        concentrations=[Concentration(Ef=1.0, carrier=CarrierConcentration(p=9.124767652322644e-06, n=9.124767651562896e-06), defects=[DefectConcentration(name='Va_O1', charges=[0], concentrations=[1.0])])])
    assert actual == expected