# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import Dict

from pymatgen import Element, IStructure

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import SingleDefectEnergy
from pydefect.input_maker.defect_entry import DefectEntry


class SingleDefectEnergyMaker:
    def __init__(self,
                 perfect: CalcResults,
                 defect: CalcResults,
                 defect_entry: DefectEntry,
                 abs_chem_pot: Dict[Element, float]):
        name = defect_entry.name
        charge = defect_entry.charge
        n_diffs = num_atom_differences(defect.structure, perfect.structure)
        energy = (defect.energy - perfect.energy
                  + reservoir_energy(n_diffs, abs_chem_pot))
        self.single_defect_energy = SingleDefectEnergy(name, charge, energy)


def reservoir_energy(diffs: Dict[Element, int],
                     abs_chem_pot: Dict[Element, float]):
    return sum([diff * abs_chem_pot[elem] for elem, diff in diffs.items()])


def num_atom_differences(structure: IStructure,
                         ref_structure: IStructure) -> Dict[Element, int]:
    c1 = structure.composition.as_dict()
    c2 = ref_structure.composition.as_dict()
    result = {}
    for k in set(c1.keys()) | set(c2.keys()):
        diff = int(c1[k] - c2[k])
        if diff:
            result[Element(k)] = diff
    return result



