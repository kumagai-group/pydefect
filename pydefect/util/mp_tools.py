# -*- coding: utf-8 -*-

from typing import List

from pymatgen import Element, MPRester

from pydefect.defaults import defaults

elements = [e.name for e in Element]


class MpQuery:
    def __init__(self,
                 element_list: List[str],
                 e_above_hull: float = defaults.e_above_hull,
                 properties: List[str] = None,
                 ):

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


