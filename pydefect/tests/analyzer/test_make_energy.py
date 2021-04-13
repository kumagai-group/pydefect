# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.energy import Energy
from pydefect.analyzer.make_energy import make_energy
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen import IStructure, Lattice


def test_make_energy(mocker):
    defect_entry = mocker.Mock(DefectEntry, autospec=True)
    defect_entry.name = "Va_Mg1"
    defect_entry.charge = -1

    calc_results = mocker.Mock(CalcResults, autospec=True)
    calc_results.structure = IStructure(Lattice.cubic(1.0), ["O"], [[0.0]*3])
    calc_results.energy = 10.0
    calc_results.electronic_conv = False

    p_calc_results = mocker.Mock(CalcResults, autospec=True)
    p_calc_results.structure = IStructure(Lattice.cubic(1.0),
                                          ["Mg", "O"], [[0.0]*3]*2)
    p_calc_results.energy = 1.0

    actual = make_energy(defect_entry, calc_results, p_calc_results)
    expected = Energy(name="Va_Mg1", charge=-1, rel_energy=10.0 - 1.0,
                      atom_io={"Mg": -1})
    assert actual == expected
