# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from typing import Dict

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import num_atom_differences
from pydefect.analyzer.energy import EnergySummary, Energy
from pydefect.corrections.abstract_correction import Correction
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen import Element


def make_energy(defect_entry: DefectEntry,
                calc_results: CalcResults,
                correction: Correction,
                perfect_calc_results: CalcResults,
                standard_energies: Dict[Element, float]) -> EnergySummary:
    atom_io = num_atom_differences(calc_results.structure,
                                   perfect_calc_results.structure)
    formation_energy = calc_results.energy - perfect_calc_results.energy
    for k, v in atom_io.items():
        formation_energy -= standard_energies[k] * v
    energy = Energy(formation_energy=formation_energy,
                    atom_io=atom_io,
                    correction_energy=correction.correction_dict)

    return EnergySummary(defect_entry.name, defect_entry.charge, energy)
