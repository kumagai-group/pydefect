# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import pydefect.tests.analyzer.test_defect_energy
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.calc_summary import SingleCalcSummary, CalcSummary
from pydefect.analyzer.make_calc_summary import make_calc_summary
from pydefect.defaults import defaults
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen import IStructure, Lattice


def test_make_calc_summary(mocker):
    defect_entry = mocker.Mock(spec=DefectEntry, autospec=True)
    defect_entry.name = "Va_O1"
    defect_entry.charge = 1
    defect_entry.full_name = "Va_O1_1"

    calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    calc_results.structure = IStructure(Lattice.cubic(1.0), ["Mg"], [[0.0]*3])
    pydefect.tests.analyzer.test_defect_energy.energy = 10.0 + defaults.abs_strange_energy - 0.1

    p_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    p_calc_results.structure = IStructure(Lattice.cubic(1.0), ["Mg", "O"], [[0.0]*3]*2)
    pydefect.tests.analyzer.test_defect_energy.energy = 10.0

    actual = make_calc_summary(defect_entry_list=[defect_entry],
                               calc_results_list=[calc_results],
                               p_calc_results=p_calc_results)

    single_summary = SingleCalcSummary(charge=1,
                                       atom_io={"O": -1},
                                       electronic_conv=calc_results.electronic_conv,
                                       ionic_conv=calc_results.ionic_conv,
                                       is_energy_strange=False)
    expected = CalcSummary({"Va_O1_1": single_summary})

    assert actual == expected


"""
TODO
-

DONE
"""
