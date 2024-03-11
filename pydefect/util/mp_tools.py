# -*- coding: utf-8 -*-

from typing import List

from pydefect.defaults import defaults
from pymatgen.core import Element

# from pymatgen.ext.matproj import MPRester, MPRestError
from mp_api.client import MPRester
from vise.util.logger import get_logger
from itertools import combinations, chain
import json
import os

elements = [e.name for e in Element]
cwd = os.path.dirname(os.path.realpath(__file__))
mp_data_file =cwd +  "/MPID_elements_10th_Jan_2024.json"

with open(mp_data_file, "r") as f:
    mp_data = json.load(f)

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

def get_compound_with_elist(element_list, mp_data_file=mp_data_file):
    with open(mp_data_file, "r") as f:
        mp_data = json.load(f)
    selected = []
    for key in mp_data.keys():
        if set(mp_data[key]).issubset(set(element_list)):
            selected.append(key)

    return selected

class MpQuery:
    def __init__(
        self,
        element_list: List[str],
        e_above_hull: float = defaults.e_above_hull,
        properties: List[str] = None,
    ):
        # API key is parsed via .pmgrc.yaml
        with MPRester() as m:
            excluded = list(set(elements) - set(element_list))
            logger.info(
                "Note that you should use the newer MPRester with package mp-api!!"
            )
            default_fields = [
                "material_id",
                "formula_pretty",
                "structure",
                "symmetry",
                "band_gap",
                "total_magnetization",
                "types_of_magnetic_species",
                "elements",
            ]
            properties = properties or default_fields
            # Adding following try and except to avoid error
            # because currently materials project is not accepting
            # excluded elements more than 15 character
            try:
                materials = m.materials.summary.search(
                    exclude_elements=excluded,
                    energy_above_hull=(-1e-5, e_above_hull),
                    fields=properties,
                )
                self.materials = materials
            except:
                # Try screening from a local database to reduce load on materials project
                try:
                    selected_ids = get_compound_with_elist(element_list,mp_data_file=mp_data_file)
                    materials = m.materials.summary.search(
                        material_ids=selected_ids, energy_above_hull=(-1e-5, e_above_hull), fields=properties
                        )
                    self.materials = materials
                except:
                    materials = m.materials.summary.search(
                        energy_above_hull=(-1e-5, e_above_hull), fields=properties
                    )
                    self.materials = exclude_rest_element(
                    docs=materials, element_list=element_list
                )


class MpEntries:
   def __init__(
       self,
       element_list: List[str],
       e_above_hull: float = defaults.e_above_hull,
       additional_properties: List[str] = None,
   ):
       excluded = list(set(elements) - set(element_list))
       criteria = {
           "elements": {"$in": element_list, "$nin": excluded},
           "e_above_hull": {"$lte": e_above_hull},
       }
       with MPRester() as m:
           try:
               materials = m.materials.summary.search(
                   exclude_elements=excluded,
                   energy_above_hull=(-1e-5, e_above_hull),
                   fields=additional_properties,
               )
               self.materials = materials
           except:
               # Try screening from a local database to reduce load on materials project
                try:
                    selected_ids = get_compound_with_elist(element_list,mp_data_file=mp_data_file)
                    materials = m.materials.summary.search(
                        material_ids=selected_ids, energy_above_hull=(-1e-5, e_above_hull), fields=additional_properties
                        )
                    self.materials = materials
                except:
                    materials = m.materials.summary.search(
                        energy_above_hull=(-1e-5, e_above_hull), fields=additional_properties
                        )
                    self.materials = exclude_rest_element(
                        docs=materials, element_list=element_list)
