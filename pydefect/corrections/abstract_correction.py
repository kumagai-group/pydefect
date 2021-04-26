# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict

from monty.json import MSONable
from monty.serialization import loadfn


class Correction(ABC, MSONable):
    def to_json_file(self, filename: str = "correction.json") -> None:
        Path(filename).write_text(self.to_json())

    @classmethod
    def from_json_file(cls, filename: str = "correction.json"):
        return loadfn(filename)
