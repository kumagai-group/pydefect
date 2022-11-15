# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from pydefect.analyzer.defect_energy import DefectEnergies, DefectEnergy
from pydefect.analyzer.defect_energy_util import u_values_from_defect_energies


def test_calc_u_value():
    defect_energies = DefectEnergies(
        atom_io={"O": -1}, charges=[0, 1, 2, 4],
        defect_energies=[
            DefectEnergy(0.0, {"correction_a": 10.0}, False),
            DefectEnergy(6.0, {"correction_a": 20.0}, False),
            DefectEnergy(4.0, {"correction_a": 30.0}, False),
            DefectEnergy(1000.0, {"correction_a": 10000.0}, True)])

    expected = {"Va_O1": {(0, 1, 2): 44.0 - 52.0}}
    assert u_values_from_defect_energies({"Va_O1": defect_energies}) == expected

    expected = {"Va_O1": {(0, 1, 2): -8.0}}
    assert u_values_from_defect_energies({"Va_O1": defect_energies},
                                         correction=False) == expected


def test_calc_u_value_with_shallow():
    defect_energies = DefectEnergies(
        atom_io={"O": -1}, charges=[0, 1, 2, 4],
        defect_energies=[
            DefectEnergy(0.0, {"correction_a": 10.0}, False),
            DefectEnergy(6.0, {"correction_a": 20.0}, False),
            DefectEnergy(4.0, {"correction_a": 30.0}, True)])

    expected = {"Va_O1": {}}
    assert u_values_from_defect_energies({"Va_O1": defect_energies}) == expected
    expected = {"Va_O1": {(0, 1, 2): 44.0 - 52.0}}
    assert u_values_from_defect_energies(
        {"Va_O1": defect_energies}, allow_shallow=True) == expected




