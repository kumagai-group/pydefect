# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass

from monty.json import MSONable

from pydefect.util.mix_in import ToJsonFileMixIn


@dataclass
class Shallow(MSONable, ToJsonFileMixIn):
    is_shallow: bool
    method: str
