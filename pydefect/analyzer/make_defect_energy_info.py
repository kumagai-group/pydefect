# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from typing import Dict, Union

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.energy import DefectEnergyInfo, DefectEnergy
from pydefect.corrections.abstract_correction import Correction
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen import Element, IStructure


def make_defect_energy_info(defect_entry: DefectEntry,
                            calc_results: CalcResults,
                            correction: Correction,
                            perfect_calc_results: CalcResults,
                            standard_energies: Dict[Element, float]
                            ) -> DefectEnergyInfo:
    atom_io = num_atom_differences(calc_results.structure,
                                   perfect_calc_results.structure)
    formation_energy = calc_results.energy - perfect_calc_results.energy
    for k, v in atom_io.items():
        formation_energy -= standard_energies[k] * v
    energy = DefectEnergy(formation_energy=formation_energy,
                          atom_io=atom_io,
                          correction_energy=correction.correction_dict)

    return DefectEnergyInfo(defect_entry.name, defect_entry.charge, energy)


def num_atom_differences(structure: IStructure,
                         ref_structure: IStructure,
                         str_key: bool = False
                         ) -> Dict[Union[Element, str], int]:
    target_composition = structure.composition.as_dict()
    reference_composition = ref_structure.composition.as_dict()
    result = {}
    for k in set(target_composition.keys()) | set(reference_composition.keys()):
        n_atom_diff = int(target_composition[k] - reference_composition[k])
        if n_atom_diff:
            result[k if str_key else Element(k)] = n_atom_diff
    return result
