# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import List, Dict
from pathlib import Path
from copy import deepcopy

from pydefect.input_maker.defect_name import DefectName
from pydefect.input_maker.supercell_info import SupercellInfo
from pydefect.database.database import oxidation_state, electronegativity


class DefectSet:
    def __init__(self, names: List[DefectName]):
        self.names = names

    def to_file(self, filename: str = "defect.in") -> None:
        Path(filename).write_text("\n".join([str(name) for name in self.names]))

    @classmethod
    def from_file(cls, filename: str = "defect.in") -> "DefectSet":
        names = []
        for line in Path(filename).read_text().split("\n"):
            names.append(DefectName.from_str(line.strip()))
        return cls(names)


class DefectSetMaker:
    def __init__(self,
                 supercell_info: SupercellInfo,
                 overwritten_oxi_states: Dict[str, int] = None,
                 dopants: List[str] = None):
        self.supercell_info = supercell_info
        self._oxidation_states = deepcopy(oxidation_state)
        if overwritten_oxi_states:
            self._oxidation_states.update(overwritten_oxi_states)

        self._dopants = dopants

        self.defect_set = self._create_defect_set()

    def _create_defect_set(self):
        result = self._create_vacancy_set()
        result.extend(self._create_substitutional_set())
        return DefectSet(result)

    def _create_vacancy_set(self):
        result = []
        for name, site in self.supercell_info.sites.items():
            ox_state = self._oxidation_states[site.element]
            for charge in charge_set(ox_state):

                result.append(DefectName("Va", name, charge))
        return result

    def _create_substitutional_set(self):
        print(electronegativity["Li"])
        print(electronegativity["H"])
        print(electronegativity["He"])
#        for dopant

        return [        DefectName.from_str("Li_H1_0"),
                        DefectName.from_str("Li_He1_-1"),
                        DefectName.from_str("Li_He1_0"),
                        DefectName.from_str("Li_He1_1"),]


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
