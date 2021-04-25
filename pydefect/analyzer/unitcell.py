# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from collections import OrderedDict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import List

import numpy as np
from monty.json import MSONable
from monty.serialization import loadfn
from ruamel.yaml import add_constructor, resolver, YAML


# courtesy of https://qiita.com/konomochi/items/f5f53ba8efa07ec5089b
add_constructor(resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                lambda loader, node: OrderedDict(loader.construct_pairs(node)))

yaml = YAML()
yaml.default_flow_style = False


@dataclass
class Unitcell(MSONable):
    system: str
    vbm: float
    cbm: float
    ele_dielectric_const: List[List[float]]
    ion_dielectric_const: List[List[float]]

    @property
    def dielectric_constant(self):
        total = (np.array(self.ele_dielectric_const)
                 + np.array(self.ion_dielectric_const))
        return total.tolist()

    def to_yaml_file(self, filename: str = "unitcell.yaml"):
        with open(Path(filename), 'w', encoding='utf-8') as f:
            yaml.dump(asdict(self), f)

    @classmethod
    def from_yaml(cls, filename: str = "unitcell.yaml"):
        return cls(**loadfn(filename))
