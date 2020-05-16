# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path

from monty.serialization import loadfn
from vise.util.logger import get_logger
logger = get_logger(__name__)

rcore = loadfn(Path(__file__).parent / "rcore.yaml")


def electronegativity(element):
    yaml_data = loadfn(Path(__file__).parent / "electronegativity.yaml")
    try:
        return yaml_data[element]
    except KeyError:
        logger.info(f"Electronegativity of {element} is set to 0.0.")
        return None


def oxidation_state(element):
    yaml_data = loadfn(Path(__file__).parent / "oxidation_state.yaml")
    try:
        return yaml_data[element]
    except KeyError:
        logger.info(f"Oxidation state of {element} is set to 0.")
        return 0
