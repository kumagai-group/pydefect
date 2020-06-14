# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import Dict, Optional

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import SingleDefectEnergy
from pydefect.corrections.abstract_correction import Correction
from pydefect.corrections.manual_correction import NoCorrection
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen import Element, IStructure


def make_single_defect_energy(perfect: CalcResults,
                              defect: CalcResults,
                              defect_entry: DefectEntry,
                              abs_chem_pot: Dict[Element, float],
                              correction: Optional[Correction] = NoCorrection()
                              ) -> SingleDefectEnergy:
    n_diffs = num_atom_differences(defect.structure, perfect.structure)
    energy = (defect.energy - perfect.energy
              + reservoir_energy(n_diffs, abs_chem_pot))
    return SingleDefectEnergy(defect_entry.name, defect_entry.charge, energy,
                              correction.correction_energy)


def reservoir_energy(diffs: Dict[Element, int],
                     abs_chem_pot: Dict[Element, float]) -> float:
    return sum([-diff * abs_chem_pot[elem] for elem, diff in diffs.items()])


def num_atom_differences(structure: IStructure,
                         ref_structure: IStructure) -> Dict[Element, int]:
    target_composition = structure.composition.as_dict()
    reference_composition = ref_structure.composition.as_dict()
    result = {}
    for k in set(target_composition.keys()) | set(reference_composition.keys()):
        n_atom_diff = int(target_composition[k] - reference_composition[k])
        if n_atom_diff:
            result[Element(k)] = n_atom_diff
    return result



