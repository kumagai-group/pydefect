# -*- coding: utf-8 -*-

from typing import List

from pydefect.defaults import defaults
from pymatgen.core import Element
#from pymatgen.ext.matproj import MPRester, MPRestError
from mp_api.client import MPRester
from vise.util.logger import get_logger
from itertools import combinations, chain

elements = [e.name for e in Element]


logger = get_logger(__name__)

def exclude_rest_element(docs, element_list):
    final_materials = []
    for i in range(len(docs)):
        elements = [e.name for e in docs[i].elements]
        if set(elements).issubset(set(element_list)):
            final_materials.append(docs[i])
        else:
            pass

    return final_materials

class MpQuery:
    def __init__(self,
                 element_list: List[str],
                 e_above_hull: float = defaults.e_above_hull,
                 properties: List[str] = None):
        # API key is parsed via .pmgrc.yaml
        with MPRester() as m:
            excluded = list(set(elements) - set(element_list))[:5]
            logger.info("Note that you should use the newer MPRester with package mp-api!!")
            default_fields = ["material_id", "formula_pretty", "structure",
                              "symmetry", "band_gap", "total_magnetization",
                              "types_of_magnetic_species","elements"]
            properties = properties or default_fields
            materials = m.materials.summary.search(
                #exclude_elements=excluded, # Currently Materials Project is not allowing excluding list longer than 15 character
                energy_above_hull=(-1e-5, e_above_hull),
                fields=properties)
            #self.materials = materials # Uncomment this when above uncomment works and comment the bottom line
            self.materials = exclude_rest_element(docs=materials, element_list=element_list)

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


