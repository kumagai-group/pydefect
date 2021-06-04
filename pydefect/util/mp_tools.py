# -*- coding: utf-8 -*-

from typing import List

from pydefect.defaults import defaults
from pymatgen.core import Element
from pymatgen.ext.matproj import MPRester

elements = [e.name for e in Element]


class MpQuery:
    def __init__(self,
                 element_list: List[str],
                 e_above_hull: float = defaults.e_above_hull,
                 properties: List[str] = None):
        default_properties = ["task_id", "full_formula", "final_energy",
                              "structure", "spacegroup", "band_gap",
                              "total_magnetization", "magnetic_type"]
        properties = properties or default_properties
        excluded = list(set(elements) - set(element_list))
        criteria = ({"elements": {"$in": element_list, "$nin": excluded},
                     "e_above_hull": {"$lte": e_above_hull}})
        # API key is parsed via .pmgrc.yaml
        with MPRester() as m:
            # Due to mp_decode=True by default, class objects are restored.
            self.materials = m.query(criteria=criteria,
                                     properties=properties)


class MpEntries:
    def __init__(self,
                 element_list: List[str],
                 e_above_hull: float = defaults.e_above_hull,
                 additional_properties: List[str] = None):
        excluded = list(set(elements) - set(element_list))
        criteria = ({"elements": {"$in": element_list, "$nin": excluded},
                     "e_above_hull": {"$lte": e_above_hull}})
        with MPRester() as m:
            self.materials = m.get_entries(
                criteria, property_data=additional_properties)


