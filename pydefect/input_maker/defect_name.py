# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from dataclasses import dataclass
import re
from typing import Union

from monty.json import MSONable


@dataclass()
class DefectName(MSONable):
    in_atom: str
    out_atom: str
    charge: int

    def __str__(self):
        return "_".join([self.in_atom, self.out_atom, str(self.charge)])

    def __hash__(self):
        return hash(str(self))

    @property
    def name_wo_charge(self):
        return "_".join([self.in_atom, self.out_atom])

    @classmethod
    def from_str(cls, string) -> "DefectName":
        name = string.split("_")
        assert len(name) == 3
        return cls(in_atom=name[0], out_atom=name[1], charge=int(name[2]))


