# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import argparse

from pymatgen.io.vasp.inputs import Poscar, Incar, Potcar
from vise.util.logger import get_logger


logger = get_logger(__name__)


def get_defect_charge_state(poscar: Poscar, potcar: Potcar, incar: Incar):
    """Get defect charge state from structure, potcar, and NELECT in INCAR. """
    nelect = incar.get("NELECT", 0)
    if poscar.site_symbols != potcar.symbols:
        raise ValueError(f"Sequence of elements in POSCAR {poscar.site_symbols}"
                         f" and that in POTCAR {potcar.symbols} is different.")
    num_elect_neutral = sum([num_atom * pot.nelectrons
                             for num_atom, pot in zip(poscar.natoms, potcar)])
    return int(num_elect_neutral - nelect)  # minus of difference of electrons


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--poscar", default="POSCAR", metavar="FILE",
                        type=Poscar.from_file, help="POSCAR file name.")
    parser.add_argument("--potcar", default="POTCAR", metavar="FILE",
                        type=Potcar.from_file, help="POTCAR file name")
    parser.add_argument("-i", "--incar", default="INCAR", metavar="FILE",
                        type=Incar.from_file, help="INCAR file name.")
    opts = parser.parse_args()

    charge_state = get_defect_charge_state(opts.poscar, opts.potcar, opts.incar)
    logger.info(f"""
Composition: {opts.poscar.structure.composition}
Charge state: {charge_state}   """)


if __name__ == "__main__":
    main()
