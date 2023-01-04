# -*- coding: utf-8 -*-

from typing import List

from pydefect.defaults import defaults
from pymatgen.core import Element
from pymatgen.ext.matproj import MPRester
from vise.util.logger import get_logger

elements = [e.name for e in Element]


logger = get_logger(__name__)


class MpQuery:
    def __init__(self,
                 element_list: List[str],
                 e_above_hull: float = defaults.e_above_hull,
                 properties: List[str] = None):
        excluded = list(set(elements) - set(element_list))
        # API key is parsed via .pmgrc.yaml
        with MPRester() as m:
            # Due to mp_decode=True by default, class objects are restored.
            try:
                default_properties = ["task_id", "full_formula", "final_energy",
                                      "structure", "spacegroup", "band_gap",
                                      "total_magnetization", "magnetic_type"]
                criteria = (
                    {"elements": {"$in": element_list, "$nin": excluded},
                     "e_above_hull": {"$lte": e_above_hull}})
                self.materials = m.query(
                    criteria=criteria,
                    properties=properties or default_properties)
            except NotImplementedError:
                logger.info("Note that you're using the newer MPRester.")
                default_fields = ["material_id", "formula_pretty", "structure",
                                  "symmetry", "band_gap", "total_magnetization",
                                  "types_of_magnetic_species"]
                self.materials = m.summary.search(
                    exclude_elements=excluded,
                    energy_above_hull=[-1e-5, e_above_hull],
                    fields=properties or default_fields)


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


