# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path
from typing import List

import numpy as np
from pydefect.analyzer.vesta.vesta_file import add_density
from pymatgen.io.vasp import Chgcar
from scipy.constants.codata import physical_constants as pc

_minor = 1e-2
default_isurfs = [0.1, 0.5, 0.8]


def make_light_chgcar(chgcar: Chgcar,
                      filename: Path,
                      max_int: int = 9,
                      vesta_file: Path = None,
                      to_vesta_file: Path = None,
                      isurfs: List[float] = None,
                      ):
    _max = np.max(chgcar.data["total"])
    total = (chgcar.data["total"] + _max * _minor) * (1 - _minor)
    new_total = np.floor((total / _max * (max_int + 1))).astype(int)
    text = chgcar.structure.to("POSCAR")
    text += f" \n {chgcar.dim[0]} {chgcar.dim[1]} {chgcar.dim[2]}\n"
    text += " ".join(new_total.flatten("F").astype(str))  # column major
    filename.write_text(text)

    if vesta_file:
        atomic_unit_to_angstrom = pc["atomic unit of length"][0] * 10 ** 10
        volume = chgcar.structure.volume / atomic_unit_to_angstrom ** 3
        isurfs = np.array(isurfs or default_isurfs) / volume * max_int
        isurfs = np.round(isurfs, 5).tolist()
        add_density(vesta_file, to_vesta_file, isurfs, filename)
