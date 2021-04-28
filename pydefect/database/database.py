# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path

from monty.serialization import loadfn
from vise.util.logger import get_logger

logger = get_logger(__name__)

electroneg_data = loadfn(Path(__file__).parent / "electronegativity.yaml")
oxi_state_data = loadfn(Path(__file__).parent / "oxidation_state.yaml")


def electronegativity(element):
    if element not in electroneg_data:
        logger.info(f"Electronegativity of {element} is set to None.")
        electroneg_data[element] = None
    return electroneg_data[element]


def oxidation_state(element):
    if element not in oxi_state_data:
        logger.info(f"Oxidation state of {element} is set to 0.")
        oxi_state_data[element] = 0
    return oxi_state_data[element]

