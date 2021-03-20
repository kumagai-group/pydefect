#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from pymatgen.io.vasp.inputs import Poscar, Incar, Potcar


def get_defect_charge_state(poscar: Poscar, potcar: Potcar, incar: Incar):
    """Get defect charge state from structure, potcar, and NELECT in INCAR. """
    nelect = incar.get("NELECT", 0)
    if poscar.site_symbols != potcar.symbols:
        raise ValueError("The sequence of elements in POTCAR and Structure "
                         "is different.")

    num_elect_neutral = sum([num_atom * pot.nelectrons
                             for num_atom, pot in zip(poscar.natoms, potcar)])
    return int(num_elect_neutral - nelect)  # minus of difference of electrons


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                        type=str,
                        help="POSCAR-type file name for the unitcell.",
                        metavar="FILE")
    parser.add_argument("--potcar", dest="potcar", default="POTCAR", type=str,
                        help="POTCAR-type file name for the unitcell.",
                        metavar="FILE")
    parser.add_argument("--incar", dest="incar", default="INCAR", type=str,
                        help="INCAR-type file name.", metavar="FILE")

    opts = parser.parse_args()

    pos = Poscar.from_file(opts.poscar)
    nions = pos.natoms
    charge = get_defect_charge_state(nions, opts.potcar, opts.incar)

    print("composition:", pos.structure.composition)
    print("defect charge:", charge)


if __name__ == "__main__":
    main()
