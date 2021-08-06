# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.calc_summary import SingleCalcSummary, CalcSummary
from pydefect.analyzer.defect_structure_info import DefectStructureInfo
from pydefect.analyzer.make_calc_summary import make_calc_summary
from pydefect.defaults import defaults
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen.core import IStructure, Lattice


def test_make_calc_summary(mocker):
    defect_entry = mocker.Mock(spec=DefectEntry, autospec=True)
    defect_entry.name = "Va_O1"
    defect_entry.charge = 1
    defect_entry.full_name = "Va_O1_1"

    calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    calc_results.structure = \
        IStructure(Lattice.cubic(1.0), ["Mg"], [[0.0]*3])
    calc_results.energy = 10.0 + defaults.abs_strange_energy - 0.1

    structure_info = mocker.Mock(spec=DefectStructureInfo, autospec=True)

    p_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    p_calc_results.structure = \
        IStructure(Lattice.cubic(1.0), ["Mg", "O"], [[0.0]*3]*2)
    p_calc_results.energy = 10.0

    actual = make_calc_summary(
        calc_set=[(calc_results, defect_entry, structure_info)],
        p_calc_results=p_calc_results)
    single_summary = SingleCalcSummary(
        charge=1,
        atom_io={"O": -1},
        electronic_conv=calc_results.electronic_conv,
        ionic_conv=calc_results.ionic_conv,
        is_energy_strange=False,
        same_config_from_init=structure_info.same_config_from_init,
        defect_type=str(structure_info.defect_type),
        symm_relation=str(structure_info.symm_relation))
    expected = CalcSummary({"Va_O1_1": single_summary})

    assert actual == expected


