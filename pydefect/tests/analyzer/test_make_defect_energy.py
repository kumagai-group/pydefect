# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import SingleDefectEnergy
from pydefect.analyzer.make_defect_energy import make_single_defect_energy, \
    num_atom_differences, reservoir_energy
from pydefect.corrections.manual_correction import ManualCorrection
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen import IStructure, Lattice, Element


def test_make_single_defect_energy(mocker):
    perfect = mocker.Mock(CalcResults, autospec=True)
    perfect.energy = 1.0
    perfect.structure = IStructure(Lattice.cubic(1), ["H"] * 2, [[0] * 3] * 2)

    defect = mocker.Mock(CalcResults, autospec=True)
    defect.energy = 3.0
    defect.structure = IStructure(Lattice.cubic(1), ["H"], [[0] * 3])

    name, charge = "Va_H1", 5
    defect_entry = mocker.Mock(DefectEntry, autospec=True)
    defect_entry.name = name
    defect_entry.charge = charge

    single_defect_energy = \
        make_single_defect_energy(perfect=perfect,
                                  defect=defect,
                                  defect_entry=defect_entry,
                                  abs_chem_pot={Element.H: 100},
                                  correction=ManualCorrection(7.0))
    energy = 3.0 - 1.0 + 100

    assert isinstance(single_defect_energy, SingleDefectEnergy)
    assert single_defect_energy == SingleDefectEnergy(name, charge, energy, 7.0)


def test_reservoir_energy():
    num_atom_diff = {Element.He: 1, Element.Li: -1}
    abs_chem_pot = {Element.He: 10.0, Element.Li: 100.0, Element.H: 1000.0}
    actual = reservoir_energy(num_atom_diff, abs_chem_pot)
    assert actual == -10.0 + 100.0


def test_num_atom_diff():
    s1 = IStructure(Lattice.cubic(1), ["H", "He"], [[0] * 3] * 2)
    s2 = IStructure(Lattice.cubic(1), ["H", "Li"], [[0] * 3] * 2)
    assert num_atom_differences(s1, s2) == {Element.He: 1, Element.Li: -1}

