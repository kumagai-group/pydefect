# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from typing import Dict

from pydefect.analyzer.band_edge_states import BandEdgeStates
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergy, DefectEnergyInfo
from pydefect.analyzer.unitcell import Unitcell
from pydefect.chem_pot_diag.chem_pot_diag import StandardEnergies
from pydefect.corrections.abstract_correction import Correction
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen import IStructure


def make_defect_energy_info(defect_entry: DefectEntry,
                            calc_results: CalcResults,
                            correction: Correction,
                            perfect_calc_results: CalcResults,
                            standard_energies: StandardEnergies,
                            unitcell: Unitcell,
                            band_edge_states: BandEdgeStates = None
                            ) -> DefectEnergyInfo:
    atom_io = num_atom_differences(calc_results.structure,
                                   perfect_calc_results.structure)

    formation_energy = calc_results.energy - perfect_calc_results.energy
    formation_energy += defect_entry.charge * unitcell.vbm
    for k, v in atom_io.items():
        formation_energy -= standard_energies[k] * v

    is_shallow = band_edge_states.is_shallow if band_edge_states else None
    energy = DefectEnergy(formation_energy=formation_energy,
                          energy_corrections=correction.correction_dict,
                          is_shallow=is_shallow)

    return DefectEnergyInfo(defect_entry.name, defect_entry.charge,
                            atom_io=atom_io, defect_energy=energy)


def num_atom_differences(structure: IStructure,
                         ref_structure: IStructure,
                         ) -> Dict[str, int]:
    target_composition = structure.composition.as_dict()
    reference_composition = ref_structure.composition.as_dict()
    result = {}
    for k in set(target_composition.keys()) | set(reference_composition.keys()):
        n_atom_diff = int(target_composition[k] - reference_composition[k])
        if n_atom_diff:
            result[k] = n_atom_diff
    return result
