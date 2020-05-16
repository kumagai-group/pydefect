# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from itertools import product
from typing import Dict, List

from pydefect.database.database import oxidation_state, electronegativity
from pydefect.input_maker.defect_name import DefectName
from pydefect.input_maker.defect_set import DefectSet
from pydefect.input_maker.supercell_info import SupercellInfo


class DefectSetMaker:
    def __init__(self,
                 supercell_info: SupercellInfo,
                 overwritten_oxi_states: Dict[str, int] = None,
                 dopants: List[str] = None,
                 ele_neg_diff: float = 2.0):
        self.supercell_info = supercell_info
        self._overwritten_oxi_states = overwritten_oxi_states or {}
        self._dopants = dopants or []
        self._ele_neg_diff = ele_neg_diff

        self.defect_set = self._create_defect_set()

    def _oxidation_state(self, element):
        try:
            return self._overwritten_oxi_states[element]
        except KeyError:
            return oxidation_state(element)

    @property
    def _host_elements(self):
        species_set = set(self.supercell_info.structure.species)
        return [str(specie) for specie in species_set]

    def _create_defect_set(self):
        result = self._create_vacancy_set()
        result.extend(self._create_substitutional_set())
        return DefectSet(result)

    def _create_vacancy_set(self):
        result = []
        for name, site in self.supercell_info.sites.items():
            oxi_state = self._oxidation_state(site.element)
            for charge in charge_set(oxi_state):
                result.append(DefectName("Va", name, charge))

        return result

    def _create_substitutional_set(self):
        result = []
        subs = self._dopants + self._host_elements
        for sub_element, (name, site) \
                in product(subs, self.supercell_info.sites.items()):
            host_element = site.element
            sub_en = electronegativity(sub_element)
            host_en = electronegativity(host_element)
            are_same_elem = sub_element == host_element
            if None in (sub_en, host_en) or are_same_elem or \
                    abs(sub_en - host_en) > self._ele_neg_diff:
                continue

            sub_oxi_state = self._oxidation_state(sub_element)
            host_oxi_state = self._oxidation_state(host_element)
            oxi_state_diff = sub_oxi_state - host_oxi_state
            for charge in charge_set(oxi_state_diff):
                result.append(DefectName(sub_element, name, charge))

        return result


def charge_set(ox_state: int) -> List[int]:
    """Set of defect charge states.

    -1 (1) is included for positive (negative) odd number.
    E.g., charge_set(3) = [-1, 0, 1, 2, 3]
          charge_set(-3) = [-3, -2, -1, 0, 1]
          charge_set(2) = [0, 1, 2]
          charge_set(-4) = [-4, -3, -2, -1, 0]

    Args:
        ox_state (int): an integer

    Returns:
        Set of candidate charges
    """
    if ox_state >= 0:
        charges = [i for i in range(ox_state + 1)]
        if ox_state % 2 == 1:
            charges.insert(0, -1)
    else:
        charges = [i for i in range(ox_state, 1)]
        if ox_state % 2 == 1:
            charges.append(1)

    return charges