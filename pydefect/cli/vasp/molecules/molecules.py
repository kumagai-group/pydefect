# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path

import yaml
from vise.util.logger import get_logger

""" 
Characteristic data is retrieved from 
McQuarrie and Simon, Phys. Chem. A molecular approach.
"""

logger = get_logger(__name__)

with open(Path(__file__).parent / "molecule_data.yaml", 'r') as f:
    MOLECULE_DATA = yaml.safe_load(f)


