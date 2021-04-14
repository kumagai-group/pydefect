# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.energy import EnergySummary, Energy
from pydefect.analyzer.make_energy_summary import make_energy_summary, \
    num_atom_differences
from pydefect.corrections.abstract_correction import Correction
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen import IStructure, Lattice, Element


def test_make_energy(mocker):
    defect_entry = mocker.Mock(DefectEntry, autospec=True)
    defect_entry.name = "Va_Mg1"
    defect_entry.charge = -1

    calc_results = mocker.Mock(CalcResults, autospec=True)
    calc_results.structure = IStructure(Lattice.cubic(1.0), ["O"], [[0.0]*3])
    calc_results.energy = 10.0
    calc_results.electronic_conv = False

    correction = mocker.Mock(Correction, autospec=True)
    correction.correction_dict = {"a": 10.0}

    p_calc_results = mocker.Mock(CalcResults, autospec=True)
    p_calc_results.structure = IStructure(Lattice.cubic(1.0),
                                          ["Mg", "O"], [[0.0]*3]*2)
    p_calc_results.energy = 1.0

    standard_energies = {Element.Mg: 10.0, Element.O: 20.0}

    actual = make_energy_summary(defect_entry, calc_results, correction, p_calc_results,
                                 standard_energies)
    expected = EnergySummary(name="Va_Mg1", charge=-1,
                             energy=Energy(formation_energy=10.0 - 1.0 + 10,
                                           atom_io={Element.Mg: -1},
                                           correction_energy={"a": 10.0}))
    assert actual == expected


def test_num_atom_diff():
    s1 = IStructure(Lattice.cubic(1), ["H", "He"], [[0] * 3] * 2)
    s2 = IStructure(Lattice.cubic(1), ["H", "Li"], [[0] * 3] * 2)
    assert num_atom_differences(s1, s2) == {Element.He: 1, Element.Li: -1}
    assert num_atom_differences(s1, s2, str_key=True) == {"He": 1, "Li": -1}