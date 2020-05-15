# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path

from monty.serialization import loadfn

electronegativity = loadfn(Path(__file__).parent / "electronegativity.yaml")
oxidation_state = loadfn(Path(__file__).parent / "oxidation_state.yaml")
rcore = loadfn(Path(__file__).parent / "rcore.yaml")
