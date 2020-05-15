# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pathlib import Path

from pydefect.input_maker.defect_name import DefectName
from monty.serialization import loadfn, dumpfn


class DefectSet:
    def __init__(self, names):
        self.names = names

    def to_file(self, filename: str = "defect.in") -> None:
        Path(filename).write_text("\n".join([str(name) for name in self.names]))

    @classmethod
    def from_file(cls, filename: str = "defect.in") -> "DefectSet":
        names = []
        for line in Path(filename).read_text().split("\n"):
            names.append(DefectName.from_str(line.strip()))
        return cls(names)


