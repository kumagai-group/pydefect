# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pymatgen import Element
from pymatgen.io.vasp import Outcar

for e in Element:
    try:
        o = Outcar(str(e) + "/OUTCAR-finish")
        print(f"{e:>2} {o.final_energy:11.8f}")
    except FileNotFoundError:
        pass
