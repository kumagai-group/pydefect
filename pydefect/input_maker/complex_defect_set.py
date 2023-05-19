# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from collections import Iterable
from dataclasses import dataclass
from typing import List, Dict, Optional, Iterator

from monty.serialization import loadfn


@dataclass
class ComplexDefect:
    """Not support interstitial related defect"""
    name: str
    replaced_sites: Dict[int, Optional[str]]
    charges: list

    # @classmethod
    # def from_dict(cls, d):
    #     return cls(name=d["name"], replaced_sites=, charges=tuple(d["charges"]))

    @property
    def str_list(self):
        return ["_".join([self.name, str(charge)]) for charge in self.charges]

    @property
    def charge_list(self):
        return [charge for charge in self.charges]


class ComplexDefectSet(Iterable):
    def __init__(self, defects: List[ComplexDefect]):
        self.defects = defects

    def __iter__(self) -> Iterator:
        yield from self.defects

    @classmethod
    def from_yaml(cls, filename: str = "defect_in.yaml") -> "ComplexDefectSet":
        d = loadfn(filename)
        names = []
        for name, dd in d.items():
            names.append(ComplexDefect(name=name, **dd))
        return cls(names)


