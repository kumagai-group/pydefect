# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from abc import ABC, abstractmethod
import json
from pathlib import Path

from monty.json import MSONable, MontyEncoder
from monty.serialization import loadfn


class Correction(ABC, MSONable):
    @property
    @abstractmethod
    def correction_energy(self):
        pass

    def to_json_file(self, filename: str) -> None:
        Path(filename).write_text(self.to_json())

    @classmethod
    def load_json(cls, filename):
        return loadfn(filename)
