# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from pydefect.analyzer.defect_energy import DefectEnergies, DefectEnergy
from pydefect.analyzer.defect_energy_util import calc_u_values


def test_calc_u_value():
    defect_energies = DefectEnergies(
        atom_io={"O": -1}, charges=[0, 1, 2, 4],
        defect_energies=[
            DefectEnergy(0.0, {"correction_a": 10.0}, True),
            DefectEnergy(6.0, {"correction_a": 20.0}, True),
            DefectEnergy(4.0, {"correction_a": 30.0}, True),
            DefectEnergy(1000.0, {"correction_a": 10000.0}, True)])

    expected = {(0, 1, 2): 44.0 - 52.0}
    assert calc_u_values(defect_energies) == expected

    expected = {(0, 1, 2): -8.0}
    assert calc_u_values(defect_energies, correction=False) == expected


