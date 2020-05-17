# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from dataclasses import dataclass
import re
from typing import Union

from monty.json import MSONable


@dataclass(frozen=True)
class Defect(MSONable):
    name: str
    charges: tuple

    @property
    def str_list(self):
        return ["_".join([self.name, str(charge)]) for charge in self.charges]

    @property
    def charge_list(self):
        return [charge for charge in self.charges]


class SimpleDefect(Defect):
    def __init__(self, in_atom, out_atom, charge_list):
        if in_atom is None:
            in_atom = "Va"
        super().__init__("_".join([in_atom, out_atom]), tuple(charge_list))

    @property
    def in_atom(self):
        result = self.name.split("_")[0]
        if result == "Va":
            return
        return result

    @property
    def out_atom(self):
        return self.name.split("_")[1]

