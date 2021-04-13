# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import num_atom_differences
from pydefect.analyzer.energy import Energy
from pydefect.input_maker.defect_entry import DefectEntry


def make_energy(defect_entry: DefectEntry,
                calc_results: CalcResults,
                perfect_calc_results: CalcResults) -> Energy:
    atom_io = num_atom_differences(calc_results.structure,
                                   perfect_calc_results.structure, True)

    return Energy(name=defect_entry.name,
                  charge=defect_entry.charge,
                  rel_energy=calc_results.energy - perfect_calc_results.energy,
                  atom_io=atom_io)
