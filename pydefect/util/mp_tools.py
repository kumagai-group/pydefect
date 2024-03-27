# -*- coding: utf-8 -*-

from typing import List

from mp_api.client import MPRester

from pydefect.defaults import defaults
from vise.util.logger import get_logger


logger = get_logger(__name__)


class MpQuery:
    def __init__(self,
                 element_list: List[str],
                 e_above_hull: float = defaults.e_above_hull,
                 properties: List[str] = None):
        # API key is parsed via .pmgrc.yaml
        with MPRester() as m:
            # Due to mp_decode=True by default, class objects are restored.
            logger.info("Note that you're using the newer MPRester.")
            default_fields = ["material_id", "formula_pretty", "structure",
                              "symmetry", "band_gap", "total_magnetization",
                              "types_of_magnetic_species"]
            properties = properties or default_fields
            self.materials = m.materials.summary.search(
                chemsys=element_list,
                energy_above_hull=(-1e-5, e_above_hull),
                fields=properties)


