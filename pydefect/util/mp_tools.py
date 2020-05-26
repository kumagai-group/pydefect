# -*- coding: utf-8 -*-

from pathlib import Path
from typing import List

from pymatgen import Element, MPRester

elements = [e.name for e in Element]
mol_dir = Path(__file__).parent / ".." / "chem_pot_diag" / "molecules"


class MpQuery:
    def __init__(self,
                 element_list: List[str],
                 properties: List[str] = None,
                 e_above_hull: float = 1e-5):

        self.element_list = element_list
        self.properties = (properties or
                           ["task_id", "full_formula", "final_energy",
                            "structure", "spacegroup", "band_gap",
                            "total_magnetization", "magnetic_type"])
        self.e_above_hull = e_above_hull

        excluded = list(set(elements) - set(element_list))
        # API key is parsed via .pmgrc.yaml

        with MPRester() as m:
            # Due to mp_decode=True by default, class objects are restored.
            self.materials = \
                m.query(criteria={"elements": {"$in": element_list,
                                               "$nin": excluded},
                                  "e_above_hull": {"$lte": e_above_hull}},
                        properties=self.properties)


