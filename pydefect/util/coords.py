# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from typing import List, Union

from vise.util.typing import Coords


def pretty_coords(coords: Union[List[float], Coords]):
    return f'({", ".join([f"{c:6.3f}" for c in coords])})'
