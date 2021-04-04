# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from math import prod
from pathlib import Path
from typing import List

import fire
import numpy as np
from pydefect.cli.vasp.make_defect_charge_info import make_spin_charges
from pymatgen.io.vasp import Chgcar, Poscar
from vise.util.logger import get_logger

logger = get_logger(__name__)


def make_spin_decomposed_parchg(parchg: Chgcar,
                                band_idx: int,
                                isosurfs: List[float] = None):
    spin_charges = make_spin_charges(parchg)

    for c, spin in zip(spin_charges, ["up", "down"]):
        if not isosurfs:
            c.write_file(f"PARCHG_{band_idx}_{spin}")
        else:
            data = np.zeros(prod(parchg.dim), dtype=int)
            normalized_vals = c.data["total"] / np.max(c.data["total"])
            for isosurf in isosurfs:
                # transpose needed as vasp is based on column measure (Fortran)
                data += (normalized_vals > isosurf).T.flatten()
            lines = [Poscar(parchg.structure).get_string()]
            lines.append(" ".join([str(d) for d in parchg.dim]))
            lines.append(" ".join(data.astype(str)))
            Path(f"PARCHG_{band_idx}_{spin}").write_text("\n".join(lines))


def make_spin_decomposed_parchg_from_files(parchg_fname: str,
                                           isosurfs: str = None):
    band_idx = int(parchg_fname.split(".")[-2])
    logger.info(f"band index {band_idx} is parsed from filename {parchg_fname}")
    if isinstance(isosurfs, float):
        isosurfs = [isosurfs]
    elif isinstance(isosurfs, str):
        isosurfs = [float(i) for i in isosurfs.split()]
    make_spin_decomposed_parchg(Chgcar.from_file(parchg_fname), band_idx, isosurfs)


if __name__ == '__main__':
    fire.Fire(make_spin_decomposed_parchg_from_files)
